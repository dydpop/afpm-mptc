function T = validate_discrete_loop()
% 闭环 Simulink 模型与数值脚本的一致性验证。
% 容差在任何仿真之前写死；不允许为了通过而放宽或改工况。
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'src','matlab'));
addpath(fileparts(mfilename('fullpath')));
outDir = fullfile(root,'data','research_v2');
if ~isfolder(outDir), mkdir(outDir); end

TOL_SPEED_RPM = 0.5;    % 转速最大绝对误差
TOL_TORQUE_RMS = 0.15;  % 电磁转矩 RMS 误差
TOL_CURRENT_RMS = 0.15; % 相电流 RMS 误差
T_END = 0.6;

logPath = fullfile(outDir,'closed_loop_validation.log');
fid = fopen(logPath,'w','n','UTF-8');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'AFPM 闭环一致性验证 | %s\n',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fprintf(fid,'MATLAB %s | Simulink license: %d\n',version,license('test','Simulink'));
fprintf(fid,'预设容差：转速 %.3f r/min，转矩 RMS %.3f N*m，电流 RMS %.3f A，时长 %.2f s\n', ...
    TOL_SPEED_RPM,TOL_TORQUE_RMS,TOL_CURRENT_RMS,T_END);

model = 'closed_loop_discrete';
if ~isfile(fullfile(root,'src','model',[model '.slx']))
    build_discrete_loop();
end
bdclose('all');
load_system(model);
set_param(model,'StopTime',sprintf('%.12g',T_END));
tStart = tic;
out = sim(model);
wall = toc(tStart);

wm_s = out.log_wm; ia_s = out.log_ia; te_s = out.log_Te;
t_s = wm_s.time; wm_s = wm_s.signals.values(:);
ia_s = ia_s.signals.values(:);
te_s = te_s.signals.values(:);

P = afpm_defaults(); P.t_end = T_END;
o = afpm_simulate('mptc',[],P);

% 时间基对齐：Simulink 首个样本在 t=0（初值），脚本首个样本在 Ts（步末）。
% 因此 sim 的第 k+1 个样本对应脚本第 k 个样本，两者逐点一一对应。
assert(abs(t_s(1))<1e-15 && abs(o.t(1)-P.Ts)<1e-12,'时间基前提不成立，需重新确认对齐方式');
wm_s = wm_s(2:end); ia_s = ia_s(2:end); te_s = te_s(2:end);
t_s = t_s(2:end);
assert(numel(t_s)==numel(o.t),'样本数不一致，不能逐点比较');
wn = o.wm; ian = o.ia; ten = o.Te;

q = t_s >= 0.05;   % 前50 ms为起动暂态，误差统计从0.05 s起（判据事先固定）
err_wm_rpm = (wm_s - wn)*60/(2*pi);
rows = { ...
    'speed_rpm', max(abs(err_wm_rpm(q))), sqrt(mean(err_wm_rpm(q).^2)), TOL_SPEED_RPM; ...
    'torque_Nm', max(abs(te_s(q)-ten(q))), sqrt(mean((te_s(q)-ten(q)).^2)), TOL_TORQUE_RMS; ...
    'current_A',  max(abs(ia_s(q)-ian(q))), sqrt(mean((ia_s(q)-ian(q)).^2)), TOL_CURRENT_RMS};
T = cell2table(rows,'VariableNames',{'signal','max_abs_err','rms_err','tolerance'});
T.passed = T.rms_err <= T.tolerance;
writetable(T,fullfile(outDir,'closed_loop_validation.csv'));

fprintf(fid,'模型运行 %.1f s；Simulink 采样点 %d；脚本采样点 %d\n',wall,numel(t_s),numel(o.t));
fprintf(fid,'统计窗口 t>=0.05 s，共 %d 点\n',nnz(q));
for k = 1:height(T)
    fprintf(fid,'%-12s max=%.3g rms=%.3g tol=%.3g -> %s\n', ...
        T.signal{k},T.max_abs_err(k),T.rms_err(k),T.tolerance(k),string(T.passed(k)));
end
fprintf(fid,'末速：Simulink %.2f r/min，脚本 %.2f r/min\n', ...
    wm_s(end)*60/(2*pi), o.wm(end)*60/(2*pi));
fprintf(fid,'结论文本只允许写"离散实现与数值脚本一致性"，不得写成接入主电路的开关级闭环。\n');

f = figure('Visible','off','Color','w','Position',[80 80 1000 620]); tiledlayout(2,1,'TileSpacing','compact');
nexttile; plot(o.t,o.wm*60/(2*pi),'k','LineWidth',1.2); hold on;
plot(t_s,wm_s*60/(2*pi),'Color',[0 .45 .74],'LineWidth',1);
ylabel('Speed (rpm)'); grid on; legend('MATLAB script','Simulink closed loop','Location','southeast');
title(sprintf('Closed-loop parity over %.1f s (max |\\Delta n| = %.2e rpm)',T_END,T.max_abs_err(1)));
nexttile; plot(o.t,o.Te,'k','LineWidth',1.2); hold on; plot(t_s,te_s,'Color',[.85 .33 .1],'LineWidth',1);
xlim([T_END-0.1 T_END]); ylabel('Torque (N m)'); xlabel('Time (s)'); grid on;
legend('MATLAB script','Simulink closed loop','Location','best');
title('Torque near the end of the run');
exportgraphics(f,fullfile(root,'figures','research_v2','closed_loop_parity.png'),'Resolution',300); close(f);

if ~all(T.passed)
    warning('闭环一致性未全部通过，详见 %s',logPath);
end
disp(T);
end
