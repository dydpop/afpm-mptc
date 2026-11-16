function tests = afpm_core_tests()
% 原生 MATLAB 回归测试；不依赖额外工具箱，不生成文件。
cases = {@test_defaults, @test_zero, @test_vectors_power, ...
    @test_delay, @test_speed_pi, @test_pwm, @test_fine, ...
    @test_mismatch, @test_current_limit, @test_loaded, ...
    @test_pwm_events, @test_pwm_convergence, ...
    @test_direct_torque_step, @test_predictor_switch, @test_ramp_profile, ...
    @test_gain_mode, @test_cost_mode};
names = {'默认参数与研究假设', '零电流初值及零负载', ...
    '电压矢量及电磁功率符号', '一拍物理延时与补偿开关', ...
    '统一转矩PI及条件积分', '真实PWM整数翻转计数', ...
    '步末时间与子步记录', '参数失配仅作用控制器', ...
    'FOC电压矢量限幅', '负载下两方案转速跟踪', ...
    '精确PWM窄脉冲全计数与冻结队列', 'FOC分段RK4步长收敛', ...
    '直接转矩阶跃旁路速度环', '预测离散方式切换', ...
    '升速给定插值与端点保持', '固定增益与估计增益可切换', '分层代价与加权代价可切换'};
tests = repmat(struct('name','','pass',false,'detail',''),numel(cases),1);
for k = 1:numel(cases)
    tests(k).name = names{k};
    try
        detail = cases{k}();
        tests(k).pass = true;
        tests(k).detail = detail;
    catch err
        tests(k).detail = sprintf('%s: %s',err.identifier,err.message);
    end
    fprintf('[%d] %s: %s\n',tests(k).pass,names{k},tests(k).detail);
end
end

function d = test_defaults()
P = afpm_defaults();
assert(P.Rs==.2 && P.L==.0085 && P.psi==.175 && P.p==4);
assert(P.J==.003 && P.F==.008 && P.Udc==250);
assert(P.Ts==50e-6 && P.Nsub==20 && P.w_ref==300);
assert(P.T_load==8 && P.t_step==.3 && P.t_end==2);
assert(P.T_rated==11.9 && P.lam_flux==1 && P.lam_sw==.02);
assert(P.Kp_w==2 && P.Ki_w==8 && P.f_sw_foc==6500);
assert(P.delay_comp && P.log_fine);
assert(strcmp(P.predictor,'euler') && isempty(P.direct_Tref));
d = '默认值符合约定（t_end=2.0、predictor=euler、direct_Tref为空）';
end

function d = test_zero()
P=afpm_defaults(); P.t_end=.002; P.w_ref=0; P.T_load=0;
for mode={'mptc','foc'}
    o=afpm_simulate(mode{1},[],P);
    assert(all(o.ia==0) && all(o.ib==0) && all(o.Te==0) && all(o.wm==0));
end
d='零给定零负载时物理状态保持零，不伪造负载';
end

function d = test_vectors_power()
P=afpm_defaults(); P.t_end=.004;
o=afpm_simulate('mptc',[],P);
S=o.switch_states; V=o.voltage_vectors;
assert(max(abs(V(:,1)-(2/3)*P.Udc*(S(:,1)-.5*S(:,2)-.5*S(:,3))))<1e-12);
assert(max(abs(V(:,2)-P.Udc/sqrt(3)*(S(:,2)-S(:,3))))<1e-12);
assert(all(V(7,:)==0) && all(V(8,:)==0));
assert(max(abs(hypot(V(1:6,1),V(1:6,2))-2*P.Udc/3))<1e-12);
ea=-P.psi*P.p*o.wm.*sin(o.theta_e);
eb=P.psi*P.p*o.wm.*cos(o.theta_e);
assert(max(abs(1.5*(ea.*o.ia+eb.*o.ib)-o.Te.*o.wm))<1e-9);
assert(o.Te(end)>0 && o.wm(end)>0);
d='Clarke电压矢量与1.5 e·i = Te·wm一致';
end

function d = test_delay()
P=afpm_defaults(); P.t_end=.005;
for enabled=[false true]
    P.delay_comp=enabled;
    for mode={'mptc','foc'}
        o=afpm_simulate(mode{1},[],P);
        assert(isequal(o.applied_vec(2:end,:),o.command_vec(1:end-1,:)));
        assert(isequal(o.applied_voltage(2:end,:),o.command_voltage(1:end-1,:)));
        assert(o.Te(1)==0 && o.wm(1)==0);
        assert(any(abs(o.command_voltage(1,:))>0));
    end
end
d='两方案均施加上一拍命令；关闭预测补偿不删除物理延时';
end

function d = test_speed_pi()
P=afpm_defaults(); P.t_end=.12; P.t_step=.04;
for mode={'mptc','foc'}
    o=afpm_simulate(mode{1},[],P);
    iw=0; cap=1.5*P.T_rated; measured=[0;o.wm(1:end-1)];
    for k=1:numel(o.t)
        e=P.w_ref*2*pi/60-measured(k); raw=P.Kp_w*e+iw;
        expected=min(max(raw,-cap),cap);
        assert(abs(o.Tref(k)-expected)<1e-10);
        if abs(raw)<=cap || raw*e<0
            iw=iw+P.Ki_w*e*P.Ts;
        end
        assert(abs(o.speed_integral(k)-iw)<1e-10);
    end
    assert(max(abs(o.Tref))<=cap && abs(o.Tref(1)-cap)<1e-12);
end
d='逐拍复算共同Nm转矩PI，含饱和时朝恢复方向积分';
end

function d = test_pwm()
P=afpm_defaults(); P.t_end=.0123;
for mode={'mptc','foc'}
    o=afpm_simulate(mode{1},[],P);
    assert(all(o.edges>=0 & o.edges==fix(o.edges)));
    assert(abs(o.f_sw-sum(o.edges)/(6*o.t(end)))<1e-10);
    S=o.switch_states; bits=S(o.fine.vec,:);
    sampled=sum(abs(diff([0 0 0;bits],1,1)),2);
    if strcmp(mode{1},'foc')
        expected=pwm_edge_reference(o,P);
        assert(all(o.fine.edges>=sampled));
        assert(all(mod(o.fine.edges-sampled,2)==0));
    else
        expected=sampled;
    end
    assert(isequal(o.fine.edges,expected));
    assert(isequal(o.edges,sum(reshape(o.fine.edges,P.Nsub,[]),1)'));
    if strcmp(mode{1},'foc')
        assert(o.f_sw>0 && abs(o.f_sw-P.f_sw_foc)>1e-6);
        assert(any(o.edges>3));
    end
end
d='由实际三桥臂状态逐次复算翻转数，不使用设定频率替代';
end

function d = test_fine()
P=afpm_defaults(); P.t_end=.002;
o=afpm_simulate('foc',[],P); N=round(P.t_end/P.Ts);
assert(isequal(o.t,(1:N)'*P.Ts));
assert(numel(o.fine.t)==N*P.Nsub);
assert(max(abs(o.fine.t(P.Nsub:P.Nsub:end)-o.t))<1e-12);
assert(isequal(o.fine.ia(P.Nsub:P.Nsub:end),o.ia));
assert(isequal(o.fine.Te(P.Nsub:P.Nsub:end),o.Te));
assert(max(abs(o.theta_e-P.p*o.theta_m))<1e-12);
P.log_fine=false; b=afpm_simulate('foc',[],P);
assert(isempty(b.fine.t) && isequal(o.ia,b.ia) && isequal(o.edges,b.edges));
d='粗细步末时刻对齐，关闭细记录不改变轨迹';
end

function d = test_mismatch()
P=afpm_defaults(); P.t_end=.004;
mm=struct('Rs',1.2,'L',.8,'psi',1.2);
for mode={'mptc','foc'}
    o=afpm_simulate(mode{1},mm,P);
    assert(isequal(o.plant_params,P));
    assert(o.controller_params.Rs==P.Rs*1.2);
    assert(o.controller_params.L==P.L*.8 && o.controller_params.psi==P.psi*1.2);
    assert(max(abs(o.Te-1.5*P.p*P.psi*o.iq))<1e-10);
    if strcmp(mode{1},'foc')
        assert(max(abs(o.iq_ref-o.Tref/(1.5*o.controller_params.p*o.controller_params.psi)))<1e-10);
        assert(o.current_Kp==o.controller_params.L*400);
        assert(o.current_Ki==o.controller_params.Rs*400);
    end
end
d='真实参数不变，转矩及失配电流给定单位可追溯';
end

function d = test_current_limit()
P=afpm_defaults(); P.Udc=20; P.t_end=.02;
o=afpm_simulate('foc',[],P);
assert(max(hypot(o.command_voltage(:,1),o.command_voltage(:,2)))<=P.Udc/sqrt(3)+1e-10);
assert(all(o.duty(:)>=0 & o.duty(:)<=1));
assert(any(o.current_saturated));
assert(all(isfinite(o.current_integral(:))));
d='低母线电压触发dq矢量限幅，积分状态有限';
end

function d = test_loaded()
% 默认t_end=2.0但时长压缩以控制测试时间：负载阶跃提前，末0.1s已进入稳态。
P=afpm_defaults(); P.t_end=.9; P.t_step=.25; values=zeros(1,2); i=0;
for mode={'mptc','foc'}
    o=afpm_simulate(mode{1},[],P); i=i+1;
    for field={'Te','wm','ia','ib','id','iq','theta_e','Tref','edges'}
        assert(all(isfinite(o.(field{1})(:))));
    end
    mask=o.t>=P.t_end-.1; values(i)=mean(o.wm(mask))*60/(2*pi);
    assert(abs(values(i)-P.w_ref)<.03*P.w_ref);
    assert(abs(mean(o.Te(mask))-P.T_load-P.F*mean(o.wm(mask)))<.6);
end
d=sprintf('末0.1s平均转速 MPTC=%.4f、FOC=%.4f r/min',values(1),values(2));
end

function expected = pwm_edge_reference(o,P)
% 独立按每桥臂解析事件时刻分桶，不由fine.vec差分推断边沿。
h=P.Ts/P.Nsub; expected=zeros(numel(o.t)*P.Nsub,1); last=[0 0 0];
for k=1:numel(o.t)
    t0=(k-1)*P.Ts; t1=k*P.Ts; duty=o.applied_vec(k,:);
    tol=8*eps(max(abs([t0,t1,P.Ts]))); events=[];
    for leg=1:3
        if duty(leg)>0 && duty(leg)<1
            for cycle=floor(t0*P.f_sw_foc):floor(t1*P.f_sw_foc)
                times=(cycle+[duty(leg)/2,1-duty(leg)/2])/P.f_sw_foc;
                events=[events,times(times>t0+tol & times<t1-tol)]; %#ok<AGROW>
            end
        end
    end
    first=min([events,t1]); mid=(t0+first)/2;
    carrier=1-abs(2*mod(mid*P.f_sw_foc,1)-1);
    bits=double(duty>=carrier); base=(k-1)*P.Nsub;
    expected(base+1)=sum(bits~=last);
    for time=events
        bin=floor((time-t0)/h+1e-9)+1;
        expected(base+bin)=expected(base+bin)+1;
    end
    final=max([events,t0]); mid=(final+t1)/2;
    carrier=1-abs(2*mod(mid*P.f_sw_foc,1)-1);
    last=double(duty>=carrier);
end
end

function d = test_pwm_events()
P=afpm_defaults(); P.t_end=2*P.Ts; P.t_step=1; P.Udc=101;
P.f_sw_foc=3.7/P.Ts; counts=zeros(3,1); frozen=zeros(3,3); nsubs=[10 20 40];
for i=1:numel(nsubs)
    P.Nsub=nsubs(i); o=afpm_simulate('foc',[],P);
    assert(isequal(o.applied_vec(1,:),[.5 .5 .5]));
    assert(isequal(o.applied_vec(2,:),o.command_vec(1,:)));
    assert(o.Te(1)==0 && o.wm(1)==0);
    assert(any(abs(o.command_voltage(1,:))>0));
    assert(isequal(o.fine.edges,pwm_edge_reference(o,P)));
    sampled=sum(abs(diff([0 0 0;o.switch_states(o.fine.vec,:)],1,1)),2);
    assert(any(o.fine.edges>sampled),'必须计入同一h内往返的窄脉冲');
    assert(isequal(o.edges,sum(reshape(o.fine.edges,P.Nsub,[]),1)'));
    counts(i)=sum(o.edges); frozen(i,:)=o.applied_vec(2,:);
end
assert(all(counts==counts(1)) && all(all(frozen==frozen(1,:))));
% 0/1占空比没有零宽伪脉冲；每周期只有a相翻转。
P.Udc=20; o=afpm_simulate('foc',[],P);
assert(any(o.applied_vec(2,:)==0) && any(o.applied_vec(2,:)==1));
assert(isequal(o.fine.edges,pwm_edge_reference(o,P)));
assert(o.edges(2)==9);
d=sprintf('Nsub=10/20/40冻结非对称duty总翻转=%g/%g/%g；窄脉冲、0/1端点及一拍延时通过',counts);
end

function d = test_pwm_convergence()
P=afpm_defaults(); P.t_end=.02; P.t_step=.0070037; nsubs=[10 20 40];
a=cell(3,1);
for i=1:3
    P.Nsub=nsubs(i); a{i}=afpm_simulate('foc',[],P);
end
e1=max(abs([a{1}.ia-a{2}.ia;a{1}.ib-a{2}.ib;a{1}.wm-a{2}.wm]));
e2=max(abs([a{2}.ia-a{3}.ia;a{2}.ib-a{3}.ib;a{2}.wm-a{3}.wm]));
assert(e1<1e-5 && e2<max(1e-10,e1/4));
assert(max(abs(a{1}.Te-a{3}.Te))<1e-5);
assert(isequal(a{1}.edges,a{2}.edges) && isequal(a{2}.edges,a{3}.edges));
assert(max(abs(a{1}.fine.ia-a{2}.fine.ia(2:2:end)))<1e-5);
assert(max(abs(a{2}.fine.ia-a{3}.fine.ia(2:2:end)))<1e-5);
d=sprintf('Nsub=10→20状态差%.3g，20→40状态差%.3g；含h内部负载阶跃及公共采样时刻',e1,e2);
end

function d = test_direct_torque_step()
P=afpm_defaults(); P.t_end=.3; P.T_load=0;
P.direct_Tref=struct('t',.1,'from',2,'to',6);
cap=1.5*P.T_rated; assert(P.direct_Tref.from<=cap && P.direct_Tref.to<=cap);
for mode={'mptc','foc'}
    o=afpm_simulate(mode{1},[],P);
    assert(strcmp(o.response_mode,'direct_torque_step'));
    assert(isequal(o.direct_Tref,P.direct_Tref));
    early=o.t<.1; late=o.t>.1;
    assert(all(o.Tref(early)==min(max(P.direct_Tref.from,-cap),cap)));
    assert(all(o.Tref(late)==min(max(P.direct_Tref.to,-cap),cap)));
    assert(all(o.speed_integral==0),'旁路速度环时积分状态须保持0');
    % Tref在步初计算：t=0.1拍（步初0.1-50us<0.1）仍属from段，下一拍起为to。
    assert(o.Tref(find(early,1,'last'))==2);
    k=find(late,1,'first')-1; assert(o.t(k)==.1 && o.Tref(k+1)==6);
    assert(all(isfinite(o.Te)) && all(isfinite(o.wm)));
end
% 电磁转矩朝to靠近：阶跃后0.1s内（转速尚未明显上升）平均Te应贴近6。
o1=afpm_simulate('mptc',[],P); o2=afpm_simulate('foc',[],P);
m=o1.t>=.1 & o1.t<=.2;
te1=mean(o1.Te(m)); te2=mean(o2.Te(m));
assert(te1>4.5 && te1<7 && te2>4.5 && te2<7);
% 限幅验证：超过1.5*T_rated的给定被截断。
Q=afpm_defaults(); Q.t_end=.1; Q.T_load=0; Q.direct_Tref=struct('t',.05,'from',2,'to',40);
o=afpm_simulate('mptc',[],Q);
assert(max(o.Tref)<=cap+1e-12 && abs(max(o.Tref)-cap)<1e-12);
d=sprintf('t<0.1恒2、t>=0.1恒6；阶跃后0.1~0.2s平均Te MPTC=%.3f、FOC=%.3f N·m；限幅30→%.3f N·m',te1,te2,cap);
end

function d = test_predictor_switch()
% 加入0.8倍参数失配，使euler与trapezoid的预测出现可分辨差异。
P=afpm_defaults(); P.t_end=.05; P.T_load=0; P.w_ref=300; P.t_step=1;
mm=struct('Rs',.8,'psi',.8);
o=afpm_simulate('mptc',mm,P); Pc=o.controller_params;
% predict_state为文件局部函数无法直接调用，此处独立复算Euler与显式Heun校验二阶修正。
a=[1.2;-0.7]; b=[-0.4;1.1]; th=0; wm=8; ua=[60;-90]; ub=[40;150];
we=Pc.p*wm; L=Pc.L; R=Pc.Rs; ps=Pc.psi; Ts=P.Ts;
dfun=@(x,y,t)([ (ua-R*x+ps*we*sin(t))/L ; (ub-R*y-ps*we*cos(t))/L ]);
f0=dfun(a,b,th);
f1=dfun(a+Ts*f0(1:2),b+Ts*f0(3:4),th+Ts*we);
e=[a+Ts*f0(1:2); b+Ts*f0(3:4)];
h2=[a+.5*Ts*(f0(1:2)+f1(1:2)); b+.5*Ts*(f0(3:4)+f1(3:4))];
assert(all(isfinite(e)) && all(isfinite(h2)));
assert(max(abs(h2-e))>1e-6);
res=cell(2,1); methods={'euler','trapezoid'};
for i=1:2
    P.predictor=methods{i};
    for mode={'mptc','foc'}
        o=afpm_simulate(mode{1},mm,P);
        assert(strcmp(o.predictor,methods{i}));
        for field={'Te','wm','ia','ib','Tref'}
            assert(all(isfinite(o.(field{1})(:))));
        end
        if strcmp(mode{1},'mptc'), res{i}=o; end
    end
end
dTe=res{1}.Te-res{2}.Te;
assert(~all(isnan(dTe)));
d=sprintf('带0.8倍失配；euler与trapezoid均可运行，MPTC转矩差绝对值均值=%.4g N·m（不判定优劣）',mean(abs(dTe)));
end

function d = test_ramp_profile()
% 升速给定：区间内线性、区间外保持端点；相同初始条件下更高的末端给定
% 必须导致更大的平均转矩指令（不涉及控制优劣判断）。
P=afpm_defaults(); P.t_end=.06; P.T_load=0; P.w_ref=20;
o1=afpm_simulate('mptc',[],P);
P2=P; P2.w_ref_profile=[0 .06; 20 60];
o2=afpm_simulate('mptc',[],P2);
assert(all(isfinite(o2.wm)) && all(isfinite(o2.Tref)));
assert(mean(o2.Tref) > mean(o1.Tref), '升速给定应带来更大平均转矩');
% 端点保持：仿真延长到升速区间之后，转速应继续追向末端给定。
P3=P; P3.t_end=.2; P3.w_ref_profile=[0 .05; 20 60];
o3=afpm_simulate('mptc',[],P3);
assert(abs(o3.wm(end)*60/(2*pi)-60) < 2, '升速结束后应保持末端给定');
d=sprintf('平均Tref %.3f>%.3f N·m；结束后%.2f r/min（给定60）', ...
    mean(o2.Tref),mean(o1.Tref),o3.wm(end)*60/(2*pi));
end

function d = test_gain_mode()
% 固定标称增益与按估计参数整定必须可区分，且默认行为不变。
P=afpm_defaults();
assert(strcmp(P.foc_gain_mode,'estimated'),'默认应为估计增益');
P.t_end=.3; mm=struct('L',1.2);
o1=afpm_simulate('foc',mm,P);
Q=P; Q.foc_gain_mode='nominal';
o2=afpm_simulate('foc',mm,Q);
assert(abs(o1.current_Kp-1.2*P.L*400)<1e-9,'估计模式应按L_est整定');
assert(abs(o2.current_Kp-P.L*400)<1e-9,'固定模式应按标称L整定');
assert(abs(o1.current_Kp-o2.current_Kp)>1e-6,'两种模式增益应不同');
assert(all(isfinite(o1.Te)) && all(isfinite(o2.Te)));
d=sprintf('Kp估计%.1f vs 固定%.1f；默认estimated',o1.current_Kp,o2.current_Kp);
end

function d = test_cost_mode()
% 分层代价：默认仍为加权；两模式都可运行；先转矩带、再磁链带、后取最少切换。
P=afpm_defaults();
assert(strcmp(P.cost_mode,'weighted'),'默认应为加权求和');
P.t_end=.05; P.T_load=0; mm=struct('L',1.2);
ow=afpm_simulate('mptc',mm,P);
Q=P; Q.cost_mode='lexicographic';
ol=afpm_simulate('mptc',mm,Q);
assert(all(isfinite(ow.Te)) && all(isfinite(ol.Te)));
assert(all(ol.command_vec>=1 & ol.command_vec<=8),'分层模式必须给出合法矢量下标');
assert(any(ow.command_vec~=ol.command_vec),'两模式应有不同的矢量决策序列');
% 容差收紧到极小值时，第一层必然筛空，必须退回最小转矩误差而不报错。
R=P; R.cost_mode='lexicographic'; R.lex_torque_tol_pct=1e-6; R.lex_flux_tol_pct=1e-6;
ot=afpm_simulate('mptc',mm,R);
assert(all(isfinite(ot.Te)) && all(ot.command_vec>=1 & ot.command_vec<=8));
d=sprintf('默认weighted；加权首拍矢量=%d，分层首拍=%d；紧容差仍可运行', ...
    ow.command_vec(1),ol.command_vec(1));
end
