function o = afpm_simulate(mode, mismatch, P)
% o = afpm_simulate('mptc'|'foc', mismatch, P)
% mismatch为控制器参数倍率结构，如struct('L',0.8,'psi',1.2)。
% 真实plant始终使用P，状态[ialpha; ibeta; wm; theta_m]，每子步RK4。
% MPTC预测采用冻结转速的一阶Euler离散近似，不等同于plant精度。
% 控制时刻k测量；pending施加于[k,k+1)，新命令从下一控制区间生效。
% applied_vec/command_vec：MPTC为1..8矢量编号；FOC为三相PWM占空比。
% FOC区间内实际离散矢量见fine.vec，vec为区间末实际矢量。
% ia是alpha(三线制相a)，ib是beta而非相b；wm单位rad/s。
% 所有t为步末时刻；Tref及command记录对应步初计算，applied对应该区间。
if nargin<1, mode='mptc'; end
if nargin<2, mismatch=[]; end
if nargin<3, P=afpm_defaults(); end
mode=validatestring(mode,{'mptc','foc'});
if ~isfield(P,'delay_comp'), P.delay_comp=true; end
if ~isfield(P,'log_fine'), P.log_fine=true; end
validateattributes(P.Ts,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(P.Nsub,{'numeric'},{'scalar','integer','positive'});
validateattributes(P.t_end,{'numeric'},{'scalar','real','finite','positive'});
for field={'Rs','L','psi','p','J','Udc','T_rated'}
    validateattributes(P.(field{1}),{'numeric'},{'scalar','real','finite','positive'});
end
Pc=P;
if ~isempty(mismatch)
    fields=fieldnames(mismatch);
    for k=1:numel(fields)
        field=fields{k};
        if ~ismember(field,{'Rs','L','psi','p','J','F','Udc'})
            error('afpm:MismatchField','失配仅支持电机及母线物理参数倍率：%s',field);
        end
        validateattributes(mismatch.(field),{'numeric'},{'scalar','real','finite','positive'});
        Pc.(field)=P.(field)*mismatch.(field);
    end
end
Ts=P.Ts; h=Ts/P.Nsub; N=round(P.t_end/Ts);
if N<1, error('afpm:Duration','仿真时长须至少包含一个控制周期。'); end
is_mptc=strcmp(mode,'mptc');
S=[1 0 0;1 1 0;0 1 0;0 1 1;0 0 1;1 0 1;0 0 0;1 1 1];
V=(2/3)*P.Udc*[S(:,1)-.5*S(:,2)-.5*S(:,3),sqrt(3)/2*(S(:,2)-S(:,3))];
Vc=V*(Pc.Udc/P.Udc);
code_to_vec=[7;5;3;4;1;6;2;8];
o.t=(1:N)'*Ts;
fields={'Te','wm','ia','ib','id','iq','theta_m','theta_e','Tref', ...
    'edges','vec','iq_ref','speed_integral'};
for k=1:numel(fields), o.(fields{k})=zeros(N,1); end
if is_mptc, nc=1; else, nc=3; end
o.applied_vec=zeros(N,nc); o.command_vec=zeros(N,nc);
o.applied_voltage=zeros(N,2); o.command_voltage=zeros(N,2);
o.current_integral=zeros(N,2); o.current_saturated=false(N,1);
o.duty=zeros(N,3);
fine_fields={'t','ia','theta_e','Te','wm','edges','vec'};
for k=1:numel(fine_fields)
    if P.log_fine, o.fine.(fine_fields{k})=zeros(N*P.Nsub,1);
    else, o.fine.(fine_fields{k})=zeros(0,1); end
end
x=zeros(4,1); iw=0; integral_dq=zeros(2,1);
pending_vec=7; pending_voltage=[0;0]; pending_duty=[.5 .5 .5];
last_bits=[0 0 0]; wref=P.w_ref*2*pi/60; Tcap=1.5*P.T_rated;
Kpi=Pc.L*400; Kii=Pc.Rs*400; Vlim=P.Udc/sqrt(3);
if ~is_mptc
    validateattributes(P.f_sw_foc,{'numeric'},{'scalar','real','finite','positive'});
    o.current_Kp=Kpi; o.current_Ki=Kii;
end
for k=1:N
    theta=Pc.p*x(4); c=cos(theta); s=sin(theta);
    % 共同速度环：输出为Nm，而非电流；允许饱和状态朝恢复方向积分。
    ew=wref-x(3); raw=P.Kp_w*ew+iw;
    Tref=min(max(raw,-Tcap),Tcap);
    if abs(raw)<=Tcap || raw*ew<0
        iw=iw+P.Ki_w*ew*Ts;
    end
    o.Tref(k)=Tref; o.speed_integral(k)=iw;
    o.iq_ref(k)=Tref/(1.5*Pc.p*Pc.psi);
    if is_mptc
        o.applied_vec(k)=pending_vec;
        applied=V(pending_vec,:)';
        xa=x(1); xb=x(2); theta_predict=Pc.p*x(4);
        if P.delay_comp
            [xa,xb,theta_predict]=predict_euler(xa,xb,theta_predict,x(3),Vc(pending_vec,1),Vc(pending_vec,2),Pc,Ts);
        end
        [xa2,xb2,theta2]=predict_euler(xa,xb,theta_predict,x(3),Vc(:,1),Vc(:,2),Pc,Ts);
        pa=Pc.L*xa2+Pc.psi*cos(theta2); pb=Pc.L*xb2+Pc.psi*sin(theta2);
        Tp=1.5*Pc.p*Pc.psi*(-xa2*sin(theta2)+xb2*cos(theta2));
        flips=sum(S~=S(pending_vec,:),2);
        cost=abs(Tref-Tp)/P.T_rated ...
            +P.lam_flux*abs(hypot(pa,pb)-Pc.psi)/Pc.psi+P.lam_sw*flips;
        [~,command_vec]=min(cost);
        command_voltage=V(command_vec,:)';
        o.command_vec(k)=command_vec;
    else
        o.applied_vec(k,:)=pending_duty; applied=pending_voltage;
        id=c*x(1)+s*x(2); iq=-s*x(1)+c*x(2);
        err=[-id;o.iq_ref(k)-iq]; we=Pc.p*x(3);
        % dq耦合补偿与永磁反电动势均用控制器参数Pc。
        uff=[-we*Pc.L*iq;we*(Pc.L*id+Pc.psi)];
        u_raw=Kpi*err+integral_dq+uff; mag=hypot(u_raw(1),u_raw(2));
        saturated=mag>Vlim;
        if saturated, u_dq=u_raw*(Vlim/mag); else, u_dq=u_raw; end
        di=Kii*Ts*err;
        % 径向限幅后的各轴条件积分：禁止继续推动该轴饱和，允许恢复。
        allow=~saturated | di.*(u_raw-u_dq)<=0;
        integral_dq=integral_dq+di.*allow;
        command_voltage=[c*u_dq(1)-s*u_dq(2);s*u_dq(1)+c*u_dq(2)];
        phase=[command_voltage(1),-.5*command_voltage(1)+sqrt(3)/2*command_voltage(2), ...
            -.5*command_voltage(1)-sqrt(3)/2*command_voltage(2)];
        zero_seq=-.5*(max(phase)+min(phase));
        command_duty=min(max(.5+(phase+zero_seq)/P.Udc,0),1);
        o.command_vec(k,:)=command_duty; o.duty(k,:)=command_duty;
        o.current_integral(k,:)=integral_dq'; o.current_saturated(k)=saturated;
    end
    o.applied_voltage(k,:)=applied'; o.command_voltage(k,:)=command_voltage';
    edges=0;
    if ~is_mptc
        % 本拍占空比冻结；交点相位为d/2、1-d/2，含区间内全部载波周期。
        t0=(k-1)*Ts; t1=k*Ts;
        duty=pending_duty(pending_duty>0 & pending_duty<1);
        phase=[duty/2,1-duty/2];
        cycles=(floor(t0*P.f_sw_foc):floor(t1*P.f_sw_foc))';
        crossing=(cycles+phase)/P.f_sw_foc;
        crossing=unique(crossing(crossing>t0 & crossing<t1));
        crossing=crossing(:);
    end
    for j=1:P.Nsub
        idx=(k-1)*P.Nsub+j; tstart=(idx-1)*h;
        if is_mptc
            actual_vec=pending_vec; bits=S(actual_vec,:);
            e=sum(bits~=last_bits); last_bits=bits; edges=edges+e;
            u=V(actual_vec,:)';
            % 若负载阶跃位于子步内部，分裂RK4步，避免跨越不连续点。
            if P.t_step>tstart && P.t_step<tstart+h
                x=rk4_step(x,u,P,0,P.t_step-tstart);
                x=rk4_step(x,u,P,P.T_load,tstart+h-P.t_step);
            else
                load_torque=double(tstart+.5*h>=P.t_step)*P.T_load;
                x=rk4_step(x,u,P,load_torque,h);
            end
        else
            tend=idx*h;
            % 只将浮点舍入范围内的交点视作记录步边界，避免零长伪片段。
            tol=8*eps(max(abs([tstart,tend,Ts])));
            internal=crossing(crossing>tstart+tol & crossing<tend-tol);
            if P.t_step>tstart && P.t_step<tend
                internal=[internal;P.t_step];
            end
            cuts=[tstart;unique(internal);tend]; e=0;
            for q=1:numel(cuts)-1
                mid=.5*(cuts(q)+cuts(q+1));
                carrier=1-abs(2*mod(mid*P.f_sw_foc,1)-1);
                bits=double(pending_duty>=carrier);
                actual_vec=code_to_vec(1+4*bits(1)+2*bits(2)+bits(3));
                flips=sum(bits~=last_bits); last_bits=bits; e=e+flips;
                u=V(actual_vec,:)';
                load_torque=double(mid>=P.t_step)*P.T_load;
                x=rk4_step(x,u,P,load_torque,cuts(q+1)-cuts(q));
            end
            edges=edges+e;
        end
        if P.log_fine
            the=P.p*x(4);
            o.fine.t(idx)=idx*h; o.fine.ia(idx)=x(1); o.fine.theta_e(idx)=the;
            o.fine.Te(idx)=1.5*P.p*P.psi*(-x(1)*sin(the)+x(2)*cos(the));
            o.fine.wm(idx)=x(3); o.fine.edges(idx)=e; o.fine.vec(idx)=actual_vec;
        end
    end
    the=P.p*x(4); ct=cos(the); st=sin(the);
    o.ia(k)=x(1); o.ib(k)=x(2); o.wm(k)=x(3);
    o.theta_m(k)=x(4); o.theta_e(k)=the;
    o.id(k)=x(1)*ct+x(2)*st; o.iq(k)=-x(1)*st+x(2)*ct;
    o.Te(k)=1.5*P.p*P.psi*o.iq(k); o.edges(k)=edges; o.vec(k)=actual_vec;
    % 直到区间结束才更新待施加命令，绝不把新命令提前施加到本拍。
    if is_mptc, pending_vec=command_vec;
    else, pending_voltage=command_voltage; pending_duty=command_duty; end
end
o.f_sw=sum(o.edges)/(6*N*Ts);
o.mode=mode; o.plant_params=P; o.controller_params=Pc;
o.switch_states=S; o.voltage_vectors=V;
o.prediction_method='冻结机械转速的一阶Euler控制预测；真实plant使用子步RK4';
o.pwm_method='零序注入两电平PWM；精确三角载波交点及负载阶跃分段RK4，固定h波形采样（非全带宽）';
o.fine_sampling='固定h步末波形；vec为末段矢量，edges累计该h内全部实际桥臂翻转；左端事件计入本步，右端事件计入下一步';
o.vector_semantics='MPTC:1..8离散编号；FOC:三相占空比；fine.vec始终为实际离散编号';
end

function [a,b,theta_next]=predict_euler(a,b,theta,wm,ua,ub,P,Ts)
we=P.p*wm;
a=a+Ts*(ua-P.Rs*a+P.psi*we*sin(theta))/P.L;
b=b+Ts*(ub-P.Rs*b-P.psi*we*cos(theta))/P.L;
theta_next=theta+Ts*we;
end

function x=rk4_step(x,u,P,load_torque,h)
k1=plant_derivative(x,u,P,load_torque);
k2=plant_derivative(x+.5*h*k1,u,P,load_torque);
k3=plant_derivative(x+.5*h*k2,u,P,load_torque);
k4=plant_derivative(x+h*k3,u,P,load_torque);
x=x+(h/6)*(k1+2*k2+2*k3+k4);
end

function dx=plant_derivative(x,u,P,load_torque)
theta=P.p*x(4); st=sin(theta); ct=cos(theta); we=P.p*x(3);
Te=1.5*P.p*P.psi*(-x(1)*st+x(2)*ct);
dx=[(u(1)-P.Rs*x(1)+P.psi*we*st)/P.L; ...
    (u(2)-P.Rs*x(2)-P.psi*we*ct)/P.L; ...
    (Te-load_torque-P.F*x(3))/P.J;x(3)];
end
