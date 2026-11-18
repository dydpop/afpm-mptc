function make_midterm_assets()
% 中期汇报专用素材生成器。
%
% 用途：只输出"两套控制器已跑通、等开关频率对齐刚完成"这一节点的证据图，
%       不输出脉动、谐波等任何结论性指标，避免提前泄露结题结果。
%
% 数据来源：data/research_v2 下已有的 protocolB_metrics.csv（100 µs / 10 kHz 对照组）
%           与 baseline_metrics.csv（50 µs / 载波校准组），只读取 f_sw_Hz 一列。
%
% 输出：figures/midterm/alignment_before_after.png

root   = fileparts(fileparts(fileparts(mfilename('fullpath'))));
outdir = fullfile(root,'figures','midterm');
datadir= fullfile(root,'data','research_v2');
if ~exist(outdir,'dir'), mkdir(outdir); end

Pb = readtable(fullfile(datadir,'protocolB_metrics.csv'),'TextType','string');
Ba = readtable(fullfile(datadir,'baseline_metrics.csv'),'TextType','string');

fsw_m_u = Pb.f_sw_Hz(Pb.mode=="mptc" & Pb.case_id=="protoB_mptc_nominal");
fsw_f_u = Pb.f_sw_Hz(Pb.mode=="foc"  & Pb.case_id=="protoB_foc_nominal");
fsw_m_a = Ba.f_sw_Hz(Ba.mode=="mptc");
fsw_f_a = Ba.f_sw_Hz(Ba.mode=="foc");
assert(numel(fsw_m_u)==1 && numel(fsw_f_u)==1 && numel(fsw_m_a)==1 && numel(fsw_f_a)==1, ...
    '开关频率取值不是单值，检查筛选条件');

fprintf('未对齐：MPTC %.1f Hz，FOC %.1f Hz，相差 %.1f 倍\n', ...
    fsw_m_u, fsw_f_u, fsw_f_u/fsw_m_u);
fprintf('对齐后：MPTC %.1f Hz，FOC %.1f Hz，相差 %.2f 倍\n', ...
    fsw_m_a, fsw_f_a, fsw_f_a/fsw_m_a);

c_m = [0.12 0.35 0.60];   % 预测控制
c_f = [0.72 0.45 0.15];   % 矢量控制
f = figure('Visible','off','Color','w','Position',[100 100 880 430]);
tiledlayout(1,2,'TileSpacing','compact','Padding','loose');

nexttile;
b = bar([fsw_m_u fsw_f_u],0.5);
b.FaceColor = 'flat'; b.CData = [c_m; c_f];
set(gca,'YScale','log','XTickLabel',{'预测控制','矢量控制'},'FontSize',10); grid on;
ylim([600 40000]);
ylabel('等效开关频率（Hz，对数轴）');
title('(a) 未对齐：控制周期 100 \mus、载波 10 kHz');
text(1,fsw_m_u*1.35,sprintf('%.0f Hz',fsw_m_u),'HorizontalAlignment','center','FontSize',10);
text(2,fsw_f_u*1.35,sprintf('%.0f Hz',fsw_f_u),'HorizontalAlignment','center','FontSize',10);
text(1.5,20000,sprintf('相差 %.1f 倍',fsw_f_u/fsw_m_u), ...
    'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');

nexttile;
b2 = bar([fsw_m_a fsw_f_a],0.5);
b2.FaceColor = 'flat'; b2.CData = [c_m; c_f];
set(gca,'XTickLabel',{'预测控制','矢量控制'},'FontSize',10); grid on;
ylim([1200 1500]); ylabel('等效开关频率（Hz）');
title('(b) 对齐后：控制周期 50 \mus、载波校准');
text(1,fsw_m_a+22,sprintf('%.0f Hz',fsw_m_a),'HorizontalAlignment','center','FontSize',10);
text(2,fsw_f_a+22,sprintf('%.0f Hz',fsw_f_a),'HorizontalAlignment','center','FontSize',10);
text(1.5,1440,sprintf('相差 %.1f%%',100*abs(fsw_f_a-fsw_m_a)/fsw_m_a), ...
    'HorizontalAlignment','center','FontSize',11,'FontWeight','bold');

exportgraphics(f,fullfile(outdir,'alignment_before_after.png'),'Resolution',300);
close(f);
fprintf('已输出 %s\n', fullfile(outdir,'alignment_before_after.png'));
end
