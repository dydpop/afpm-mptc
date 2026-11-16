function testResult=afpm_metric_tests()
% 真实 MATLAB 合成信号测试，不调用 Python、不写出附加文件。
names=["非整数时间窗纯正弦";"五次谐波10百分比";"直流不改变THD"; ...
    "基波近零无效";"完整周期不足无效";"Nyquist不足不虚构频谱"; ...
    "fine优先及控制带宽标注";"翻转计数手算";"短暂穿越不得稳定"; ...
    "末端不足20ms不得稳定";"稳定后恢复可测";"非单调相位无效"; ...
    "瞬时速度不稳无效";"显式窗口不被覆盖";"宽带残差不等于THD"; ...
    "40次谐波幅值";"负向旋转";"边沿窗口不对齐拒绝估计"; ...
    "强高次弱基波仍无效";"无负载阶跃不报告恢复";"窗口后扰动不污染响应"; ...
    "窗口外停转重复相位不影响稳态THD"; ...
    "稳态信号门禁通过";"上升拖尾判未收敛"; ...
    "稳态偏差判偏差";"纹波无趋势仍判稳态";"全频段THD与低次并列"];
passed=false(numel(names),1); detail=strings(numel(names),1);
for k=1:numel(names)
    try
        detail(k)=run_case(k); passed(k)=true;
    catch err
        detail(k)=string(err.identifier)+": "+string(err.message);
    end
    fprintf('[%d] %s — %s\n',passed(k),names(k),detail(k));
end
testResult=table(names,passed,detail,'VariableNames',{'name','passed','detail'});
end

function d=run_case(k)
[o,P]=fixture(); d="ok";
switch k
    case 22
        o.theta_e(o.t<.1)=0;
        o.ia=sin(o.theta_e);
        m=afpm_metrics(o,P);
        assert(m.thd_valid && m.thd_pct<1e-4);
        d="startup duplicate phase excluded from valid steady window";
    case 1
        m=afpm_metrics(o,P);
        assert(m.thd_valid && m.thd_pct<1e-4);
        assert(abs(m.electrical_cycles-round(m.electrical_cycles))>.1);
        assert(m.thd_used_cycles==5 && m.thd_used_highest_harmonic==40);
        assert(m.thd_window_s(1)>=P.metric_window_s(1) && m.thd_window_s(2)<=P.metric_window_s(2));
        d=sprintf('THD=%.9g%%；相位窗口%.6f..%.6fs',m.thd_pct,m.thd_window_s);
    case 2
        o.ia=sin(o.theta_e)+.1*sin(5*o.theta_e); m=afpm_metrics(o,P);
        assert(m.thd_valid && abs(m.thd_pct-10)<.003);
        d=sprintf('THD=%.8f%%',m.thd_pct);
    case 3
        o.ia=sin(o.theta_e)+.1*sin(5*o.theta_e); a=afpm_metrics(o,P);
        o.ia=o.ia+13.7; b=afpm_metrics(o,P);
        assert(abs(a.thd_pct-b.thd_pct)<1e-8 && b.thd_valid);
    case 4
        o.ia=2+1e-12*sin(o.theta_e); m=afpm_metrics(o,P);
        assert(~m.thd_valid && isnan(m.thd_pct) && contains(m.thd_reason,'基波'));
    case 5
        P.metric_window_s=[1.05 1.2]; m=afpm_metrics(o,P);
        assert(~m.thd_valid && isnan(m.thd_pct) && contains(m.thd_reason,'5个'));
    case 6
        o=thin(o,50); m=afpm_metrics(o,P);
        assert(~m.thd_valid && isnan(m.thd_pct) && contains(m.thd_reason,'Nyquist'));
    case 7
        fine=o; o=thin(o,50); o.fine=fine;
        a=afpm_metrics(o,P); assert(a.thd_valid && strcmp(a.spectrum_source,'fine'));
        assert(strcmp(a.switching_source,'fine'));
        o.fine=struct('t',[]); b=afpm_metrics(o,P);
        assert(~b.thd_valid && strcmp(b.spectrum_source,'control'));
        assert(contains(b.bandwidth_note,'控制步'));
        assert(a.source_nyquist_Hz>b.source_nyquist_Hz);
    case 8
        o.edges=mod((1:numel(o.t))',4);
        m=afpm_metrics(o,P); expected=sum(o.edges(o.t>.8 & o.t<=1.2))/(6*.4);
        assert(abs(m.f_sw_Hz-expected)<1e-8);
        o.edges(:)=1; m=afpm_metrics(o,P);
        assert(abs(m.f_sw_Hz-1/(6*P.Ts))<1e-6);
        d=sprintf('f_sw=%.8fHz，手算同值',m.f_sw_Hz);
    case 9
        rpm=P.w_ref*ones(size(o.t)); rpm(o.t>=.3)=.9*P.w_ref;
        rpm(o.t>=.4 & o.t<.41)=P.w_ref; o.wm=rpm*2*pi/60;
        o.Te(:)=8; o.Te(o.t>=.3)=9.5; o.Te(o.t>=.4 & o.t<.41)=8;
        % 最终窗口均值8Nm，但末端交替出带，不能把首次穿越判稳。
        post=o.t>=.8; o.Te(post)=8+1.5*sin(2*pi*10*o.t(post));
        m=afpm_metrics(o,P);
        assert(isnan(m.load_recovery_ms) && ~m.load_stable_at_end);
        assert(isnan(m.torque_settle_ms) && ~m.torque_stable_at_end);
    case 10
        o.wm(o.t>=.3 & o.t<1.19)=.9*P.w_ref*2*pi/60;
        m=afpm_metrics(o,P);
        assert(isnan(m.load_recovery_ms) && ~m.load_stable_at_end);
    case 11
        o.wm(o.t>=.3 & o.t<.6)=.9*P.w_ref*2*pi/60;
        o.Te(o.t>=.3 & o.t<.55)=10;
        m=afpm_metrics(o,P);
        assert(m.load_stable_at_end && abs(m.load_recovery_ms-300)<.1);
        assert(m.torque_stable_at_end && m.torque_settle_ms>=250 && m.torque_settle_ms<252);
        assert(abs(m.load_speed_dip_rpm-.1*P.w_ref)<1e-9);
        assert(m.torque_filter_length_s>=.001-1e-12);
        assert(~m.direct_torque_step_measured);
        d=sprintf('负载恢复%.4fms，转矩稳定%.4fms',m.load_recovery_ms,m.torque_settle_ms);
    case 12
        ix=find(o.t>1,1); o.theta_e(ix)=o.theta_e(ix-1)-.1;
        m=afpm_metrics(o,P); assert(~m.thd_valid && contains(m.thd_reason,'单调'));
    case 13
        wob=1+.03*sin(2*pi*10*o.t); o.wm=o.wm.*wob;
        o.theta_e=.37+cumtrapz(o.t,P.p*o.wm); o.ia=sin(o.theta_e);
        m=afpm_metrics(o,P);
        % 3%纹波在窗口均值上无趋势，稳态门禁按趋势判定为稳态；
        % THD 仍须无效。稳态门禁优先于频谱判据，故原因可能指向稳态或2%速度变化；
        % 只要不是静默返回一个数即可。
        assert(~m.thd_valid && isnan(m.thd_pct));
        assert(contains(m.thd_reason,'2%') || contains(m.thd_reason,'准稳态'));
    case 14
        P.metric_window_s=[.7 1.1]; m=afpm_metrics(o,P);
        assert(isequal(m.metric_window_s,[.7 1.1]));
        assert(m.thd_valid && m.thd_window_s(2)<=1.1);
        P.metric_window_s=.4; m=afpm_metrics(o,P);
        assert(max(abs(m.metric_window_s-[.8 1.2]))<1e-12);
    case 15
        o.ia=sin(o.theta_e)+.2*sin(70*o.theta_e); m=afpm_metrics(o,P);
        assert(m.thd_valid && m.thd_pct<.001);
        assert(abs(m.wideband_distortion_pct-20)<.05);
        d=sprintf('THD=%.8f%%，宽带残差=%.6f%%',m.thd_pct,m.wideband_distortion_pct);
    case 16
        o.ia=sin(o.theta_e)+.1*sin(40*o.theta_e); m=afpm_metrics(o,P);
        assert(m.thd_valid && abs(m.thd_pct-10)<.1);
        d=sprintf('40次THD=%.6f%%；Nyquist=%.1fHz',m.thd_pct,m.source_nyquist_Hz);
    case 17
        P.w_ref=-P.w_ref; o.wm=-o.wm; o.theta_e=-o.theta_e;
        o.ia=sin(o.theta_e)+.1*sin(5*o.theta_e); m=afpm_metrics(o,P);
        assert(m.thd_valid && abs(m.thd_pct-10)<.003);
    case 18
        P.metric_window_s=[.800001 1.2]; m=afpm_metrics(o,P);
        assert(isnan(m.f_sw_Hz) && contains(m.f_sw_reason,'对齐'));
    case 19
        o.ia=1e-12*sin(o.theta_e)+sin(5*o.theta_e); m=afpm_metrics(o,P);
        assert(~m.thd_valid && isnan(m.thd_pct) && contains(m.thd_reason,'基波'));
    case 20
        P.T_load=0; m=afpm_metrics(o,P);
        assert(isnan(m.load_recovery_ms) && isnan(m.torque_settle_ms));
        assert(~m.direct_torque_step_measured);
    case 21
        P.metric_window_s=[.7 1.1];
        o.wm(o.t>1.15)=0; o.Te(o.t>1.15)=0; m=afpm_metrics(o,P);
        assert(m.load_stable_at_end && m.torque_stable_at_end);
        assert(m.response_window_s(2)==1.1);
    case 23
        o.wm=o.wm.*(1+.0005*sin(2*pi*10*o.t)); m=afpm_metrics(o,P);
        assert(m.steady_valid && m.ripple_valid);
        assert(strcmp(m.steady_reason,'ok') && strcmp(m.ripple_reason,'ok'));
        assert(m.thd_valid && isfinite(m.thd_pct));
        d=sprintf('趋势%.3g%%，偏差%.3g%%，THD=%.3g%%', ...
            m.speed_trend_pct,m.speed_offset_pct,m.thd_pct);
    case 24
        % 窗口内从294线性升到299 r/min，模拟速度环慢极点拖尾。
        % 必须清掉 fine 记录，否则门禁仍读到原稳态细采样，测不到拖尾。
        o.fine=struct('t',[]);
        ramp=2*pi/60*interp1([0 P.t_end],[294 299],o.t,'linear');
        o.wm=ramp;
        o.theta_e=.37+2*pi*P.p/2/pi*cumtrapz(o.t,ramp);
        m=afpm_metrics(o,P);
        assert(~m.steady_valid);
        assert(contains(m.steady_reason,'趋势未收敛'));
        assert(m.speed_trend_pct>m.steady_trend_pct_threshold);
        assert(~m.thd_valid && isnan(m.thd_pct) && contains(m.thd_reason,'准稳态'));
        assert(~m.ripple_valid && isfinite(m.ripple_pp_Nm));
        d=sprintf('趋势%.4g%%，ripple保留%.4gNm',m.speed_trend_pct,m.ripple_pp_Nm);
    case 25
        o.fine=struct('t',[]);
        o.wm=(297.4*2*pi/60)*ones(size(o.t)); o.theta_e=.37+2*pi*(P.p*297.4/60)*o.t;
        m=afpm_metrics(o,P);
        assert(~m.steady_valid && contains(m.steady_reason,'稳态偏差'));
        assert(m.speed_trend_pct<m.steady_trend_pct_threshold);
        assert(~isnan(m.speed_offset_pct) && m.speed_offset_pct>m.steady_offset_pct_threshold);
        assert(isnan(m.thd_pct) && ~m.thd_valid && ~m.ripple_valid);
        d=sprintf('偏差%.4g%%，趋势%.3g%%',m.speed_offset_pct,m.speed_trend_pct);
    case 26
        % 300±3 r/min 正弦纹波：极差达2%，但无趋势。
        wob=2*pi/60*(300+3*sin(2*pi*8*o.t));
        o.wm=wob;
        o.theta_e=.37+2*pi*P.p/2/pi*cumtrapz(o.t,wob);
        P.w_ref=300; m=afpm_metrics(o,P);
        assert(m.steady_valid);
        assert(m.speed_trend_pct<m.steady_trend_pct_threshold);
        assert(m.ripple_valid && m.thd_valid);
        d=sprintf('极差%.3g%%，趋势%.4g%% → 稳态', ...
            100*(max(wob)-min(wob))/mean(wob),m.speed_trend_pct);
    case 27
        % 全频段THD：必须与低次口径并列，且不改变原 thd_pct 的定义。
        m=afpm_metrics(o,P);
        assert(m.thd_valid && isfinite(m.thd_pct));
        assert(isfinite(m.thd_full_pct) && m.thd_full_pct>=m.thd_pct-1e-9);
        assert(numel(m.harmonic_peak_A)==40 && numel(m.harmonic_peak_full_A)==200);
        assert(isequal(m.thd_full_harmonics,2:200));
        % 纯正弦(仅基波)时两者都应接近0，且全频段不小于低次。
        o2=o; o2.theta_e=o.theta_e;
        o2.ia=sin(o.theta_e);
        m2=afpm_metrics(o2,P);
        assert(isfinite(m2.thd_full_pct) && m2.thd_full_pct<0.5);
        d=sprintf('2-40次%.3f%%，2-200次%.3f%%，开关频段%.3f%%', ...
            m.thd_pct,m.thd_full_pct,m.switching_band_pct);
end
end

function [o,P]=fixture()
P=struct('Ts',50e-6,'w_ref',317,'p',4,'t_step',.3,'t_end',1.2, ...
    'T_rated',11.9,'T_load',8,'metric_window_s',[.8 1.2]);
t=(1:round(P.t_end/P.Ts))'*P.Ts;
o=struct('t',t,'theta_e',.37+2*pi*(P.p*P.w_ref/60)*t, ...
    'wm',P.w_ref*2*pi/60*ones(size(t)),'Te',8*ones(size(t)), ...
    'edges',zeros(size(t)));
o.ia=sin(o.theta_e);
end

function o=thin(o,n)
fields=fieldnames(o);
for k=1:numel(fields), o.(fields{k})=o.(fields{k})(n:n:end); end
end
