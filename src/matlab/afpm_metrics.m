function m = afpm_metrics(o,P)
% 可审计指标。wm 为 rad/s，w_ref 为 rpm；t 是区间末时刻。
% metric_window_s 接受 [起点 终点] 或末端窗口长度，默认末 0.4 s。
% THD：整数电周期内、去直流后的 2..40 次相位谐波 RMS/基波 RMS。
% 对源样点的分段线性相位波形精确积分，不将插值 FFT 当成新增带宽。
required={'Ts','w_ref','p','t_step','t_end','T_rated','T_load'};
for k=1:numel(required)
    assert(isfield(P,required{k}) && isscalar(P.(required{k})) && ...
        isfinite(P.(required{k})), 'afpm:Parameter','参数缺失或无效：%s',required{k});
end
assert(P.Ts>0 && P.p>0 && P.t_end>0 && P.T_rated>0,'afpm:Parameter','参数范围无效');
w=[max(0,P.t_end-.4),P.t_end];
if isfield(P,'metric_window_s') && ~isempty(P.metric_window_s)
    v=P.metric_window_s;
    if isscalar(v), w=[P.t_end-v,P.t_end]; else, w=v(:)'; end
end
assert(numel(w)==2 && all(isfinite(w)) && w(1)>=0 && w(2)>w(1) ...
    && w(2)<=P.t_end+1e-10,'afpm:Window','指标窗口必须在仿真范围内');
m=struct('speed_mean_rpm',NaN,'speed_error_pct',NaN,'torque_mean_Nm',NaN, ...
    'ripple_pp_Nm',NaN,'ripple_pct',NaN,'ripple_valid',false, ...
    'ripple_reason','尚未计算','thd_pct',NaN,'thd_valid',false, ...
    'thd_reason','尚未计算','wideband_distortion_pct',NaN,'f_sw_Hz',NaN, ...
    'load_recovery_ms',NaN,'load_speed_dip_rpm',NaN,'torque_settle_ms',NaN, ...
    'metric_window_s',w,'electrical_cycles',0);
m.metric_window_source='默认末0.4s';
if isfield(P,'metric_window_s') && ~isempty(P.metric_window_s), m.metric_window_source='P.metric_window_s'; end
m.ripple_definition='峰峰值 / P.T_rated * 100';
m.speed_error_definition='(平均转速 - 给定转速) / abs(给定转速) * 100';
[tw,wm,sw]=channel(o,'wm'); [tt,Te,st]=channel(o,'Te');
m.speed_source=sw; m.torque_source=st;
[xw,yw]=clip(tw,wm,w); [xt,yt]=clip(tt,Te,w);
m.speed_mean_rpm=trapz(xw,yw)/diff(w)*60/(2*pi);
if P.w_ref~=0, m.speed_error_pct=100*(m.speed_mean_rpm-P.w_ref)/abs(P.w_ref); end
% 稳态门禁：只做判定，不改已有公式。趋势项衡量窗口内漂移，
% 避免用极差把转速纹波误判为未收敛。
m.steady_trend_pct_threshold=0.2; m.steady_offset_pct_threshold=0.5;
m.steady_definition=['窗口内转速最小二乘趋势相对变化<=0.2% ' ...
    '且 |均值-给定|/|给定|<=0.5%'];
m.speed_trend_pct=NaN; m.speed_offset_pct=NaN; m.speed_trend_raw_pct=NaN;
m.steady_valid=false; m.steady_reason='尚未判定';
rpmw=yw*60/(2*pi);
if diff(w)<=0 || ~all(isfinite(rpmw))
    m.steady_reason='窗口长度或转速数据无效';
else
    % 纯最小二乘在纹波周期数不足时会把纹波吸进斜率。先去掉OLS直线，
    % 再从两侧残差均值差反推纹波贡献并扣除，保留真实拖尾斜率。
    xr=xw-xw(1); L=xr(end);
    c=polyfit(xr,rpmw,1);
    res=rpmw-polyval(c,xr);
    ripple_gain=0;
    if L>0
        third=L/3;
        i1=xr<=third; i2=xr>=L-third;
        if any(i1) && any(i2)
            ripple_gain=abs(mean(res(i1))-mean(res(i2)))/L;
        end
    end
    slope=max(0,abs(c(1))-ripple_gain);
    mu=mean(rpmw);
    if abs(mu)>1e-12
        m.speed_trend_raw_pct=100*abs(c(1))*L/abs(mu);
        m.speed_trend_pct=100*slope*L/abs(mu);
    end
    if P.w_ref~=0, m.speed_offset_pct=100*abs(mu-P.w_ref)/abs(P.w_ref); end
    trend_bad=~isfinite(m.speed_trend_pct) || m.speed_trend_pct>m.steady_trend_pct_threshold;
    offset_bad=isnan(m.speed_offset_pct) || m.speed_offset_pct>m.steady_offset_pct_threshold;
    if trend_bad && offset_bad
        m.steady_reason='趋势未收敛且存在稳态偏差';
    elseif trend_bad
        m.steady_reason='趋势未收敛';
    elseif offset_bad
        m.steady_reason='存在稳态偏差';
    else
        m.steady_valid=true; m.steady_reason='ok';
    end
end
m.torque_mean_Nm=trapz(xt,yt)/diff(w);
m.ripple_pp_Nm=max(yt)-min(yt); m.ripple_pct=100*m.ripple_pp_Nm/P.T_rated;
% 峰峰值在未收敛窗口可能含拖尾漂移，只标记有效性，不清零数值。
if m.steady_valid
    m.ripple_valid=true; m.ripple_reason='ok';
else
    m.ripple_valid=false; m.ripple_reason='稳态窗口内未收敛；峰峰值可能含拖尾漂移';
end
[te,edges,se]=channel(o,'edges');
m.switching_source=se; m.switching_window_s=w;
m.switching_definition='三桥臂翻转数 / (6 * 窗口秒数)';
m.f_sw_reason='ok'; m.switching_edges=NaN;
tol=1e-9*max(1,P.t_end);
if any(edges<0 | edges~=fix(edges))
    m.f_sw_reason='翻转计数不是非负整数';
elseif ~(w(1)==0 || any(abs(te-w(1))<=tol)) || ~any(abs(te-w(2))<=tol)
    m.f_sw_reason='窗口边界未与计数区间对齐，不能拆分整数边沿';
else
    m.switching_edges=sum(edges(te>w(1)+tol & te<=w(2)+tol));
    m.f_sw_Hz=m.switching_edges/(6*diff(w));
end
m=spectrum(m,o,P,w);
% 准稳态门禁联动：未通过稳态判定时不能保留旧 THD 值。
% 若spectrum已给出更具体的拒绝原因（如2%判据），保留原因为不可覆盖的诊断。
if ~m.steady_valid
    m.thd_valid=false; m.thd_pct=NaN;
    if contains(m.thd_reason,'准稳态') || strcmp(m.thd_reason,'ok') ...
            || strcmp(m.thd_reason,'尚未计算')
        m.thd_reason='稳态窗口内未满足准稳态（见steady_reason）';
    end
end
m.response_kind='负载扰动恢复；不是直接转矩给定阶跃';
m.direct_torque_step_measured=false; m.response_hold_s=.020;
m.response_window_s=[P.t_step,min([w(2),tw(end),tt(end)])];
m.load_stable_at_end=false; m.torque_stable_at_end=false;
m.load_recovery_reason='未观察到有效负载阶跃'; m.torque_settle_reason=m.load_recovery_reason;
m.load_speed_dip_definition='正向给定减最小转速；反向按旋转方向折算';
m.torque_filter_samples=max([1,ceil(.001/median(diff(tt))),ceil(P.Ts/median(diff(tt)))]);
m.torque_filter_length_s=m.torque_filter_samples*median(diff(tt));
m.torque_filter_definition='因果尾随移动平均，至少1ms且不少于一个控制周期';
m.torque_target_Nm=m.torque_mean_Nm; m.torque_tolerance_Nm=.05*abs(m.torque_mean_Nm);
if P.t_step<0 || P.t_step>=m.response_window_s(2) || P.T_load==0, return; end
post=tw>=P.t_step & tw<=m.response_window_s(2);
if P.w_ref==0
    m.load_recovery_reason='零转速给定无法定义相对2%带宽';
elseif any(post)
    rpm=wm*60/(2*pi);
    m.load_speed_dip_rpm=max(0,abs(P.w_ref)-min(sign(P.w_ref)*rpm(post)));
    [m.load_recovery_ms,m.load_stable_at_end,m.load_recovery_reason]= ...
        suffix(tw(post),abs(rpm(post)-P.w_ref)<=.02*abs(P.w_ref),P.t_step,.020);
end
if abs(m.torque_mean_Nm)<1e-9
    m.torque_settle_reason='最终转矩均值近零，无法定义相对5%带宽'; return;
end
if max(abs(diff(tt)-median(diff(tt))))>1e-6*median(diff(tt))
    m.torque_settle_reason='转矩时间采样不均匀，样点移动平均不代表固定时长'; return;
end
n=m.torque_filter_samples; filt=movmean(Te,[n-1,0]);
post=tt>=P.t_step & tt<=m.response_window_s(2) & (1:numel(tt))'>=n;
[m.torque_settle_ms,m.torque_stable_at_end,m.torque_settle_reason]= ...
    suffix(tt(post),abs(filt(post)-m.torque_mean_Nm)<=m.torque_tolerance_Nm,P.t_step,.020);
end

function [t,y,source]=channel(o,field)
source='control'; d=o;
if isfield(o,'fine') && isstruct(o.fine) && isfield(o.fine,'t') ...
        && ~isempty(o.fine.t) && isfield(o.fine,field) && ~isempty(o.fine.(field))
    d=o.fine; source='fine';
end
assert(isfield(d,'t') && isfield(d,field),'afpm:Data','缺少数据字段：%s',field);
t=d.t(:); y=d.(field)(:);
assert(numel(t)>=2 && numel(t)==numel(y) && all(isfinite(t)) ...
    && all(isfinite(y)) && all(diff(t)>0),'afpm:Data','时间或数据无效：%s',field);
end

function [x,y]=clip(t,v,w)
assert(t(1)<=w(1)+1e-10 && t(end)>=w(2)-1e-10,'afpm:Coverage', ...
    '数据不能覆盖固定窗口；不外推缺失的初始状态');
inside=t>w(1) & t<w(2);
x=[w(1);t(inside);w(2)];
y=[interp1(t,v,w(1),'linear');v(inside);interp1(t,v,w(2),'linear')];
end

function m=spectrum(m,o,P,w)
m.thd_harmonics=2:40; m.thd_highest_harmonic=40;
m.thd_used_highest_harmonic=NaN; m.thd_required_cycles=5; m.thd_used_cycles=0;
m.thd_reconstruction_note='分段线性重建在接近Nyquist处会衰减；Nyquist通过不是幅值精度保证';
m.thd_method='整数电周期同步相位积分；源波形分段线性解析积分，无插值频谱外推';
m.thd_window_s=[NaN NaN]; m.thd_available_cycles=0;
m.source_nyquist_Hz=NaN; m.source_sample_rate_Hz=NaN;
m.source_max_resolvable_harmonic=NaN; m.electrical_frequency_Hz=NaN;
m.speed_variation_pct=NaN; m.phase_speed_variation_pct=NaN;
m.fundamental_peak_A=NaN; m.harmonic_peak_A=nan(1,40);
m.wideband_reason='THD先决条件未满足';
m.wideband_definition='原始采样带宽内、相位加权的去直流去基波残差RMS/基波RMS；不等于THD';
m.spectrum_source='control'; d=o;
fields={'t','ia','theta_e','wm'};
if isfield(o,'fine') && isstruct(o.fine) && all(cellfun(@(f)isfield(o.fine,f) && ~isempty(o.fine.(f)),fields))
    d=o.fine; m.spectrum_source='fine';
end
m.bandwidth_note='原始fine采样带宽；无抗混叠滤波记录时不保证模拟信号无混叠';
if strcmp(m.spectrum_source,'control')
    m.bandwidth_note='仅控制步采样；不代表开关宽带频谱，且不能排除采样前混叠';
end
if ~all(cellfun(@(f)isfield(d,f) && ~isempty(d.(f)),fields))
    m.thd_reason='缺少电流、相位或速度数据'; return;
end
t=d.t(:); ia=d.ia(:); theta=d.theta_e(:); wm=d.wm(:);
if numel(t)<3 || any([numel(ia),numel(theta),numel(wm)]~=numel(t)) ...
        || any(~isfinite([t;ia;theta;wm])) || any(diff(t)<=0)
    m.thd_reason='频谱源数据非有限、长度不一致或时间不递增'; return;
end
use=find(t>=w(1) & t<=w(2));
if numel(use)<3, m.thd_reason='窗口内采样不足'; return; end
% 加入跨越窗口边界的原始样点，以原始最大采样间隔核查带宽。
i0=max(1,use(1)-1); i1=min(numel(t),use(end)+1);
dt=diff(t(i0:i1)); dq=diff(theta(i0:i1));
m.source_nyquist_Hz=1/(2*max(dt)); m.source_sample_rate_Hz=1/median(dt);
signq=sign(median(dq));
if signq==0 || any(signq*dq<=0)
    m.thd_reason='电角度不严格单调（含停转、反转或未展开相位）'; return;
end
f=abs(dq)./(2*pi*dt);
m.electrical_frequency_Hz=mean(f);
m.source_max_resolvable_harmonic=floor(min(m.source_nyquist_Hz/max(f),pi/max(abs(dq))));
q=signq*theta;
[~,qc]=clip(t,q,w); [~,wc]=clip(t,wm,w);
m.electrical_cycles=diff(qc([1 end]))/(2*pi);
if abs(mean(wc))<1e-9
    m.thd_reason='窗口平均机械转速近零'; return;
end
m.speed_variation_pct=100*(max(wc)-min(wc))/abs(mean(wc));
m.phase_speed_variation_pct=100*(max(f)-min(f))/mean(f);
if m.speed_variation_pct>=2 || m.phase_speed_variation_pct>=2
    m.thd_reason='瞬时转速或相位速度峰峰变化达到2%，不满足准稳态'; return;
end
% 同时检查相位导数和机械转速，防止记录相位与电机参数不一致。
fwm=abs(interp1(t,wm,(t(i0:i1-1)+t(i0+1:i1))/2))*P.p/(2*pi);
if max(abs(f-fwm))/mean(f)>.02
    m.thd_reason='相位导数与p*wm不一致超过2%'; return;
end
first=ceil(qc(1)/(2*pi)-1e-10); last=floor(qc(end)/(2*pi)+1e-10);
m.thd_available_cycles=max(0,last-first);
if m.thd_available_cycles<5
    m.thd_reason='固定窗口不足5个完整电周期'; return;
end
if m.source_nyquist_Hz<=40*max(f) || max(abs(dq))*40>=pi
    m.thd_reason='原始时间/相位采样Nyquist不足以辨识40次谐波'; return;
end
qa=(last-5)*2*pi; qb=last*2*pi;
% Interpolate only inside the validated monotonic window. Startup can have
% repeated phase at standstill and must not enter the interpolation domain.
qw=q(i0:i1); iw=ia(i0:i1); tw=t(i0:i1);
inside=qw>qa & qw<qb;
x=[qa;qw(inside);qb]; y=[interp1(qw,iw,qa);iw(inside);interp1(qw,iw,qb)];
m.thd_window_s=[interp1(qw,tw,qa),interp1(qw,tw,qb)];
m.thd_used_cycles=5; m.thd_used_highest_harmonic=40;
% 相位同步的整数周期域中直接积分，避免1024点降采样引入额外混叠。
x=x-qa; L=x(end); dc=trapz(x,y)/L; y=y-dc;
dx=diff(x); mid=(x(1:end-1)+x(2:end))/2;
yavg=(y(1:end-1)+y(2:end))/2; slope=diff(y)./dx;
coef=zeros(1,40);
for h=1:40
    z=h*dx/2; s=sin(z)./z; g=(z.*cos(z)-sin(z))./z.^2;
    small=abs(z)<1e-3; zz=z(small);
    g(small)=-zz/3+zz.^3/30-zz.^5/840;
    coef(h)=2/L*sum(exp(-1i*h*mid).*(yavg.*dx.*s+1i*slope.*dx.^2/2.*g));
end
m.harmonic_peak_A=abs(coef); m.fundamental_peak_A=abs(coef(1));
ac=sqrt(trapz(x,y.^2)/L);
if m.fundamental_peak_A<max(1e-9,1e-6*ac)
    m.thd_reason='基波近零或相对交流RMS过弱'; return;
end
m.thd_pct=100*sqrt(sum(abs(coef(2:40)).^2))/abs(coef(1));
m.thd_valid=true; m.thd_reason='ok';
residual=y-real(coef(1)*exp(1i*x));
m.wideband_distortion_pct=100*sqrt(trapz(x,residual.^2)/L)/(abs(coef(1))/sqrt(2));
m.wideband_reason='ok；仅代表记录采样带宽内的残差';
end

function [ms,stable,reason]=suffix(t,inside,tstep,hold)
ms=NaN; stable=false; reason='末端未连续保持稳定20ms';
if isempty(t) || ~inside(end), return; end
lastbad=find(~inside,1,'last');
if isempty(lastbad), k=1; else, k=lastbad+1; end
% 必须持续到记录末端，短暂进入容差带不能作为恢复时间。
if t(end)-t(k)<hold-1e-10, return; end
stable=true; ms=max(0,t(k)-tstep)*1000; reason='ok';
end
