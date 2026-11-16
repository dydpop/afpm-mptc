function export_report_tables()
% 把结果 CSV 汇总成一个报告用工作簿，便于后续写结题报告/PPT 时直接取数。
% 只读 CSV、不参与实验指纹，改这里不需要重跑仿真。
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
out = fullfile(root,'data','research_v2');
xls = fullfile(out,'report_tables.xlsx');
if isfile(xls), delete(xls); end
B = readtable(fullfile(out,'baseline_metrics.csv'),'TextType','string');
S = readtable(fullfile(out,'sensitivity_summary.csv'),'TextType','string');
W = readtable(fullfile(out,'sweep_metrics.csv'),'TextType','string');
C = readtable(fullfile(out,'convergence_summary.csv'),'TextType','string');
R = readtable(fullfile(out,'robustness_metrics.csv'),'TextType','string');
D = readtable(fullfile(out,'direct_step_metrics.csv'),'TextType','string');
P = readtable(fullfile(out,'ramp_metrics.csv'),'TextType','string');
A = string(splitlines(fileread(fullfile(out,'acceptance_summary.txt'))));
A = A(strlength(strtrim(A))>0);

base = B(:,{'mode','speed_mean_rpm','speed_error_pct','ripple_pp_Nm','thd_pct', ...
    'wideband_distortion_pct','f_sw_Hz','load_recovery_ms','steady_valid', ...
    'thd_valid','foc_Kp_if_used','foc_Ki_if_used'});
sens = S(:,{'param','mode','nominal_ripple_Nm','max_ripple_change_pct', ...
    'max_ripple_change_Nm','max_thd_change_percentage_points','min_f_sw_Hz', ...
    'max_f_sw_Hz','valid_off_nominal_thd_count'});
sweep = W(:,{'case_id','mode','param','factor','speed_mean_rpm','ripple_pp_Nm', ...
    'thd_pct','f_sw_Hz','load_recovery_ms','id_control_peak_abs_A','controller_psi', ...
    'steady_valid'});
conv = C(:,{'mode','metric','N20','N40','relative_change_pct','screening_passed'});
abl = R(R.family~="integration_step",:);
abl = abl(:,{'case_id','family','mode','param','factor','speed_mean_rpm', ...
    'ripple_pp_Nm','thd_pct','f_sw_Hz','load_recovery_ms','steady_valid'});
step = D(:,{'mode','rise_ms','overshoot_pct','final_mean_Nm','final_err_pct', ...
    'sustained_20ms','f_sw_Hz'});
ramp = P(:,{'mode','initial_speed_rpm','max_track_err_rpm','rms_track_err_rpm', ...
    'peak_torque_Nm','peak_torque_pct_rated','peak_current_A','f_sw_Hz','final_speed_rpm'});
acc = table(A,'VariableNames',{'acceptance'});

writetable(base,xls,'Sheet','baseline');
writetable(sens,xls,'Sheet','sensitivity');
writetable(sweep,xls,'Sheet','sweep_30rows');
writetable(conv,xls,'Sheet','convergence');
writetable(abl,xls,'Sheet','ablations');
writetable(step,xls,'Sheet','direct_step');
writetable(ramp,xls,'Sheet','speed_ramp');
writetable(acc,xls,'Sheet','acceptance');
fprintf('已导出 %s（8 个工作表）\n',xls);
end
