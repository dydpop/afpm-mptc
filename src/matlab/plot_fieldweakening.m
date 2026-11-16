function plot_fieldweakening()
% 弱磁延伸图：三组母线/转速组合下两方案能否维持给定。
% 只读已保存的轨迹与本机 CSV，不参与实验指纹，改图不需要重跑仿真。
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(root,'data','research_v2');
figdir = fullfile(root,'figures','research_v2');
if ~isfolder(figdir), mkdir(figdir); end
D = readtable(fullfile(out,'fieldweakening_metrics.csv'),'TextType','string');
modes = ["mptc","foc"]; cols = [0 0.447 0.741; 0.85 0.325 0.098];
cases = unique(D(:,{'Udc_V','w_ref_rpm'}),'stable');
f = figure('Visible','off','Color','w','Position',[60 60 1200 640]);
tiledlayout(1,numel(D.Udc_V)/2,'TileSpacing','compact');
for i = 1:height(cases)
    nexttile; hold on; leg = {};
    for m = 1:2
        row = D(D.mode==modes(m) & D.Udc_V==cases.Udc_V(i) & D.w_ref_rpm==cases.w_ref_rpm(i),:);
        if isempty(row), continue; end
        file = fullfile(out,char(row.case_id(1))+".mat");
        if ~isfile(file), continue; end
        z = load(file); o = z.o;
        plot(o.t,o.wm*60/(2*pi),'Color',cols(m,:),'LineWidth',1.1);
        leg{end+1} = sprintf('%s (%.0f r/min, %+.1f%%)',upper(modes(m)), ...
            row.final_speed_rpm(1),row.speed_err_pct(1)); %#ok<AGROW>
    end
    yline(cases.w_ref_rpm(i),'k--','LineWidth',1);
    grid on; xlabel('Time (s)'); ylabel('Speed (rpm)');
    title(sprintf('%d V bus, reference %d r/min',cases.Udc_V(i),cases.w_ref_rpm(i)));
    if ~isempty(leg), legend([leg,{'Reference'}],'Location','southeast','FontSize',7); end
end
sgtitle('Voltage-limited operation: neither scheme has a field-weakening module','FontSize',11);
exportgraphics(f,fullfile(figdir,'field_weakening.png'),'Resolution',300); close(f);
fprintf('已生成 figures/research_v2/field_weakening.png\n');
end
