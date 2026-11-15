function mptc_sim(task)
if nargin < 1, task = 'demo'; end
% 表贴式 PMSM + MPTC / FOC 数值仿真（MATLAB，不依赖 Simulink 图形界面）
%
% 用途：
%   1) 搭 Simulink 前先验证 MPTC 算法逻辑
%   2) 参数扫描快速出数据（比 Simulink 快 1—2 个数量级）
%   3) 给建模同学提供"标准答案"对照
%
% 相对上一版（Python）修正的三处错误：
%   - 反电动势符号：psi_a = L*ia + psi*cos(th)，故 L*dia/dt = u - R*ia + psi*p*wm*sin(th)
%   - 额定转矩：0.5 HP / 300 r/min 折算应为 11.9 N·m（此前误设 2.0）
%   - 转速外环加抗积分饱和：饱和时停止积分

clc;
P = default_params();
if strcmp(task,'sweep')
    run_sweep_impl(P);
    return;
end

fprintf('=== 单工况验证：MPTC，无参数偏差 ===\n');
o = simulate('mptc', [], P);
m = calc_metrics(o, P);
fprintf('  末速 %.1f r/min（给定 %d）\n', o.wm(end)*60/(2*pi), P.w_ref);
fprintf('  稳态转矩 %.2f N·m（负载 %.1f）\n', mean(o.Te(end-3000:end)), P.T_load);
fprintf('  响应时间 %.1f ms | 脉动 %.2f N·m (%.1f%%) | THD %.1f%% | f_sw %.0f Hz\n', ...
    m.t_rise*1e3, m.ripple, m.ripple/P.T_rated*100, m.thd, m.f_sw);

fprintf('\n=== FOC 对照 ===\n');
of = simulate('foc', [], P);
mf = calc_metrics(of, P);
fprintf('  末速 %.1f r/min | 稳态转矩 %.2f | 脉动 %.2f | THD %.1f%% | f_sw %.0f Hz\n', ...
    of.wm(end)*60/(2*pi), mean(of.Te(end-3000:end)), mf.ripple, mf.thd, mf.f_sw);

fprintf('\n=== 采样周期敏感性（MPTC）===\n');
fprintf('  %8s %10s %10s %10s %10s\n','Ts(us)','末速','脉动','脉动%','f_sw');
for Ts = [100e-6 50e-6]
    P2 = P; P2.Ts = Ts;
    o2 = simulate('mptc', [], P2);
    m2 = calc_metrics(o2, P2);
    fprintf('  %8.1f %10.1f %10.2f %9.1f%% %10.0f\n', ...
        Ts*1e6, o2.wm(end)*60/(2*pi), m2.ripple, m2.ripple/P.T_rated*100, m2.f_sw);
end

plot_result(o, P, 'MPTC 基准工况');
end


% ==================== 参数 ====================
function P = default_params()
P.Rs   = 0.2;        % 定子电阻 Ω
P.L    = 8.5e-3;     % d/q 电感 H（表贴式 Ld = Lq）
P.psi  = 0.175;      % 永磁磁链 Wb
P.p    = 4;          % 极对数
P.J    = 0.003;      % 转动惯量 kg·m²
P.F    = 0.008;      % 阻尼 N·m·s
P.Udc  = 250;        % 母线电压 V
P.Ts   = 50e-6;      % 控制器采样周期 s
P.Nsub = 10;         % 每个 Ts 内的积分子步
P.w_ref  = 300;      % 转速给定 r/min
P.T_load = 8.0;      % 负载转矩 N·m
P.t_step = 0.3;      % 负载突加时刻 s
P.t_end  = 0.8;      % 仿真时长 s
P.T_rated = 11.9;    % 额定转矩（0.5HP/300rpm 折算）
P.lam_flux = 0.05;   % 磁链项权重
P.lam_sw   = 0.02;   % 开关次数项权重
P.Kp_w = 2.0;        % 转速环 P
P.Ki_w = 8.0;        % 转速环 I
P.f_sw_foc = 6500;   % FOC 载波频率 Hz（与 MPTC 实测等效频率对齐）
end


% ==================== 电压矢量 ====================
function [ua, ub, S] = voltage_vectors(Udc)
S = [1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 1 0 1; 0 0 0; 1 1 1];
sa = S(:,1); sb = S(:,2); sc = S(:,3);
ua = (2/3)*Udc*(sa - 0.5*sb - 0.5*sc);
ub = (2/3)*Udc*(sqrt(3)/2)*(sb - sc);
end


% ==================== 磁链与转矩 ====================
function [ps_abs, Te, pa, pb] = flux_torque(ia, ib, thm, P)
the = P.p * thm;
pa = P.L*ia + P.psi*cos(the);
pb = P.L*ib + P.psi*sin(the);
Te = 1.5 * P.p * (pa*ib - pb*ia);
ps_abs = hypot(pa, pb);
end


% ==================== 预测（多子步）====================
function [ia, ib, thm] = predict(ia, ib, ua, ub, wm, thm, P, h, nsub)
hs = h / nsub;
for k = 1:nsub
    the = P.p * thm;
    ea =  P.psi * P.p * wm * sin(the);
    eb = -P.psi * P.p * wm * cos(the);
    ia = ia + hs*(ua - P.Rs*ia + ea)/P.L;
    ib = ib + hs*(ub - P.Rs*ib + eb)/P.L;
    thm = thm + hs*wm;
end
end


% ==================== 主仿真 ====================
function o = simulate(mode, mismatch, P)
Udc = P.Udc; Ts = P.Ts; Nsub = P.Nsub; h = Ts/Nsub;
N = round(P.t_end / Ts);
[ua_v, ub_v, S] = voltage_vectors(Udc);

% 控制器内部参数（可与实际不同，用于模拟参数失配）
Pc = P;
if ~isempty(mismatch)
    fn = fieldnames(mismatch);
    for i = 1:numel(fn)
        Pc.(fn{i}) = P.(fn{i}) * mismatch.(fn{i});
    end
end

x = [0; 0; 0; 0];             % [ia; ib; wm; thm]
w_ref_rad = P.w_ref * 2*pi/60;

prev_vec = 7;                 % 上一拍实际施加的矢量索引（初值零矢量 000）
int_w = 0; int_d = 0; int_q = 0;
sw_count = 0;
prev_bits = S(7,:);

o.t = zeros(N,1); o.Te = zeros(N,1); o.wm = zeros(N,1);
o.id = zeros(N,1); o.iq = zeros(N,1); o.ia = zeros(N,1); o.ib = zeros(N,1);
o.vec = zeros(N,1); o.Tref = zeros(N,1);

for k = 1:N
    t = (k-1)*Ts;
    ia = x(1); ib = x(2); wm = x(3); thm = x(4);
    the = P.p * thm;
    T_now = (t >= P.t_step) * P.T_load;

    if strcmp(mode,'mptc')
        % --- 转速外环（抗积分饱和）---
        w_err = w_ref_rad - wm;
        T_max = 1.5*P.T_rated;
        T_raw = P.Kp_w*w_err + int_w;
        T_ref = min(max(T_raw, -T_max), T_max);
        if T_raw == T_ref
            int_w = min(max(int_w + P.Ki_w*w_err*Ts, -T_max), T_max);
        end

        % --- 延时补偿：用上一拍实际矢量推进一拍 ---
        [ia1, ib1, thm1] = predict(ia, ib, ua_v(prev_vec), ub_v(prev_vec), wm, thm, Pc, Ts, Nsub);

        % --- 8 矢量枚举 ---
        cost = zeros(8,1); Te_pred = zeros(8,1);
        for i = 1:8
            [ia2, ib2, thm2] = predict(ia1, ib1, ua_v(i), ub_v(i), wm, thm1, Pc, Ts, Nsub);
            [ps2, Te2] = flux_torque(ia2, ib2, thm2, Pc);
            Te_pred(i) = Te2;
            n_sw = sum(S(i,:) ~= prev_bits);
            cost(i) = abs(T_ref - Te2)/P.T_rated ...
                    + P.lam_flux*abs(Pc.psi - ps2)/Pc.psi ...
                    + P.lam_sw*n_sw;
        end
        [~, sel] = min(cost);
        ua_app = ua_v(sel); ub_app = ub_v(sel);
        sw_count = sw_count + sum(S(sel,:) ~= prev_bits);
        prev_bits = S(sel,:);
        prev_vec = sel;

    else   % FOC
        w_err = w_ref_rad - wm;
        iq_ref = min(max(2.0*w_err + int_w, -2*P.T_rated), 2*P.T_rated);
        int_w = min(max(int_w + 8.0*w_err*Ts, -2*P.T_rated), 2*P.T_rated);
        id_ref = 0;
        id_ =  ia*cos(the) + ib*sin(the);
        iq_ = -ia*sin(the) + ib*cos(the);
        % 电流环增益按控制器"认为"的 L、Rs 整定；参数失配时整定值偏离实际
        Kp_i = Pc.L * 400;
        Ki_i = Pc.Rs * 400;
        % 反电动势前馈用控制器"认为"的 psi
        uq_ff = Pc.psi * Pc.p * wm;
        ud = Kp_i*(id_ref - id_) + int_d;
        uq = Kp_i*(iq_ref - iq_) + int_q + uq_ff;
        int_d = min(max(int_d + Ki_i*(id_ref-id_)*Ts, -Udc), Udc);
        int_q = min(max(int_q + Ki_i*(iq_ref-iq_)*Ts, -Udc), Udc);
        ua_app = cos(the)*ud - sin(the)*uq;
        ub_app = sin(the)*ud + cos(the)*uq;
        mag = hypot(ua_app, ub_app);
        lim = Udc/sqrt(3);
        if mag > lim
            ua_app = ua_app*lim/mag; ub_app = ub_app*lim/mag;
        end
        sw_count = sw_count + 6*(P.f_sw_foc*Ts);   % 按目标载波频率折算开关次数
        T_ref = iq_ref;
    end

    % --- 系统积分 ---
    for j = 1:Nsub
        [dx, ~] = deriv(x, ua_app, ub_app, P, T_now);
        x = x + h*dx;
    end

    % --- 记录 ---
    [~, Te, ~, ~] = flux_torque(x(1), x(2), x(4), P);
    o.t(k) = t; o.Te(k) = Te; o.wm(k) = x(3);
    o.ia(k) = x(1); o.ib(k) = x(2);
    o.id(k) =  x(1)*cos(P.p*x(4)) + x(2)*sin(P.p*x(4));
    o.iq(k) = -x(1)*sin(P.p*x(4)) + x(2)*cos(P.p*x(4));
    o.vec(k) = prev_vec; o.Tref(k) = T_ref;
end
o.f_sw = sw_count / (6 * P.t_end);
end


function [dx, Te] = deriv(x, ua, ub, P, T_load)
ia = x(1); ib = x(2); wm = x(3); thm = x(4);
the = P.p * thm;
ea =  P.psi * P.p * wm * sin(the);
eb = -P.psi * P.p * wm * cos(the);
dia = (ua - P.Rs*ia + ea)/P.L;
dib = (ub - P.Rs*ib + eb)/P.L;
[~, Te] = flux_torque(ia, ib, thm, P);
dwm = (Te - T_load - P.F*wm)/P.J;
dx = [dia; dib; dwm; wm];
end


% ==================== 指标 ====================
function m = calc_metrics(o, P)
t = o.t; Te = o.Te; Ts = P.Ts;
i0 = find(t >= P.t_step, 1);
T_ss = mean(Te(end-3000:end));          % 以稳态值为目标，而非额定值
band = 0.05 * abs(T_ss);
idx = find(abs(Te(i0:end) - T_ss) <= band, 1);
if isempty(idx), m.t_rise = NaN; else, m.t_rise = idx*Ts; end
ss = find(t >= t(end) - 0.25, 1);
seg = Te(ss:end);
m.ripple = max(seg) - min(seg);
sig = o.ia(ss:end) - mean(o.ia(ss:end));
n = numel(sig);
if n > 256
    A = abs(fft(sig)); A = A(1:floor(n/2));
    f = (0:numel(A)-1)' * (1/Ts) / n;
    f_e = P.p * mean(o.wm(ss:end)) / (2*pi);      % 电频率，已知
    [~, k0] = min(abs(f - f_e));                   % 直接锁到基波
    if k0 >= 2 && A(k0) > 0
        m.thd = sqrt(sum(A(2:end).^2) - A(k0)^2) / A(k0) * 100;
    else
        m.thd = NaN;
    end
else
    m.thd = NaN;
end
m.f_sw = o.f_sw;
end


function plot_result(o, P, ttl)
figure('Color','w','Position',[100 100 900 620]);
subplot(3,1,1);
plot(o.t, o.wm*60/(2*pi), 'b', 'LineWidth', 1.2); hold on;
yline(P.w_ref, 'r--'); grid on; ylabel('转速 (r/min)');
title(ttl); legend('实际','给定','Location','southeast');
subplot(3,1,2);
plot(o.t, o.Te, 'b', 'LineWidth', 0.8); hold on;
plot(o.t, o.Tref, 'r--', 'LineWidth', 1.0);
yline(P.T_rated, 'k:'); grid on; ylabel('转矩 (N·m)');
legend('电磁转矩','转矩给定','额定','Location','southeast');
subplot(3,1,3);
plot(o.t, o.ia, 'b'); hold on; plot(o.t, o.ib, 'r'); grid on;
ylabel('相电流 (A)'); xlabel('时间 (s)'); legend('i_a','i_b');
outdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'figures');
if ~exist(outdir,'dir'), mkdir(outdir); end
exportgraphics(gcf, fullfile(outdir,'mptc_baseline.png'), 'Resolution', 150);
fprintf('  图已保存: figures/mptc_baseline.png\n');
end


% ==================== 参数扫描 ====================
function run_sweep_impl(P)
fprintf('=== 参数敏感性扫描 ===\n');
params = {'Rs','L','psi'};
names  = {'定子电阻 R_s','电感 L','永磁磁链 \psi_f'};
offs   = [0.8 1.0 1.2];            % 相对标称值
modes  = {'mptc','foc'};

R = table();
fprintf('  %-8s %-6s %-6s %10s %10s %10s %10s\n', ...
    '参数','方案','偏移','响应(ms)','脉动(N·m)','THD(%)','f_sw(Hz)');
for mi = 1:numel(modes)
    for pi_ = 1:numel(params)
        for oi = 1:numel(offs)
            mm = [];
            if offs(oi) ~= 1
                mm = struct(params{pi_}, offs(oi));
            end
            o = simulate(modes{mi}, mm, P);
            m = calc_metrics(o, P);
            R = [R; table({params{pi_}}, {modes{mi}}, offs(oi), ...
                m.t_rise*1e3, m.ripple, m.thd, m.f_sw, ...
                'VariableNames', {'param','mode','offset','t_rise_ms','ripple','thd','f_sw'})];
            fprintf('  %-8s %-6s %-6.2f %10.2f %10.2f %10.1f %10.0f\n', ...
                params{pi_}, modes{mi}, offs(oi), m.t_rise*1e3, m.ripple, m.thd, m.f_sw);
        end
    end
end

outdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
if ~exist(outdir,'dir'), mkdir(outdir); end
writetable(R, fullfile(outdir,'sweep_results.csv'));
fprintf('\n  数据已存: data/sweep_results.csv（%d 行）\n', height(R));

% 绘图：每个参数一张图，分两个子图（转矩脉动 / 阶跃响应时间），两方案分色
figdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'figures');
if ~exist(figdir,'dir'), mkdir(figdir); end
figure('Color','w','Position',[80 60 1200 320]);
for pi_ = 1:numel(params)
    subplot(1,3,pi_);
    selM = strcmp(R.param, params{pi_}) & strcmp(R.mode,'mptc');
    selF = strcmp(R.param, params{pi_}) & strcmp(R.mode,'foc');
    dM = R(selM,:); dF = R(selF,:);
    yyaxis left;
    plot(dM.offset, dM.ripple, '-o', 'LineWidth',1.6, 'Color',[0.85 0.33 0.10]); hold on;
    plot(dF.offset, dF.ripple, '-s', 'LineWidth',1.6, 'Color',[0.00 0.45 0.74]);
    ylabel('转矩脉动 (N\cdotm)');
    yyaxis right;
    plot(dM.offset, dM.t_rise_ms, '--o', 'LineWidth',1.1, 'Color',[0.85 0.33 0.10]);
    plot(dF.offset, dF.t_rise_ms, '--s', 'LineWidth',1.1, 'Color',[0.00 0.45 0.74]);
    ylabel('阶跃响应时间 (ms)');
    xlabel('参数倍率（1.0 = 标称）'); title(names{pi_}); grid on;
    xlim([0.75 1.25]); xticks(offs); xticklabels({'0.8','1.0','1.2'});
    if pi_ == 1
        legend({'MPTC 脉动','FOC 脉动','MPTC 响应','FOC 响应'}, ...
            'Location','northoutside','Orientation','horizontal','FontSize',8);
    end
end
sgtitle('参数失配下的性能退化（实线：转矩脉动，虚线：阶跃响应时间）','FontSize',11);
exportgraphics(gcf, fullfile(figdir,'sweep_sensitivity.png'), 'Resolution', 160);
fprintf('  图已存: figures/sweep_sensitivity.png\n');
end
