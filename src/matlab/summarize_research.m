function S = summarize_research()
% Deterministic experiment summaries; no p-values or fictitious replicates.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'data','research_v2');
B=readtable(fullfile(out,'baseline_metrics.csv'),'TextType','string');
R=readtable(fullfile(out,'sweep_metrics.csv'),'TextType','string');
C=readtable(fullfile(out,'robustness_metrics.csv'),'TextType','string');
T=readtable(fullfile(out,'core_tests.csv'),'TextType','string');
U=readtable(fullfile(out,'metric_tests.csv'),'TextType','string');
V=readtable(fullfile(out,'simulink_validation.csv'),'TextType','string');
batch=load(fullfile(out,'manifest.mat'),'manifest');
assert(strcmp(batch.manifest.status,'complete'),'Batch is not complete.');
for tab={B,R,C}
    D=tab{1};
    assert(all(D.run_id==string(batch.manifest.run_id)));
    assert(all(D.source_fingerprint==string(batch.manifest.source_fingerprint)));
    for i=1:height(D)
        file=fullfile(out,char(D.case_id(i)+".mat"));
        assert(isfile(file),'Missing raw trajectory: %s',file);
        z=load(file,'P','m');
        assert(strcmp(z.P.run_id,batch.manifest.run_id) && strcmp(z.m.source_fingerprint,batch.manifest.source_fingerprint),'Stale raw file: %s',file);
    end
end
assert(height(R)==30 && numel(unique(R.case_id))==30);
assert(all(T.passed) && all(U.passed) && all(V.passed));
assert(all(isfinite(R.speed_mean_rpm)) && all(isfinite(R.ripple_pp_Nm)) && all(isfinite(R.f_sw_Hz)));
assert(all(abs(R.speed_error_pct)<2),'A sweep case fails the mean-speed gate.');
% Complete means measured or explicitly invalid, never silently replaced by zero.
assert(all(R.thd_valid | (isnan(R.thd_pct) & strlength(R.thd_reason)>0)));
rows=cell(0,12); pars=["Rs","L","psi"]; modes=["mptc","foc"];
for ip=1:3
    for im=1:2
        D=R(R.param==pars(ip)&R.mode==modes(im),:);
        assert(isequal(sort(D.factor),[.8;.9;1;1.1;1.2]));
        nominal=D(D.factor==1,:); off=D(D.factor~=1,:);
        rp=(off.ripple_pp_Nm/nominal.ripple_pp_Nm-1)*100;
        td=(off.thd_pct/nominal.thd_pct-1)*100;
        rows(end+1,:)={pars(ip),modes(im),nominal.ripple_pp_Nm,max(rp),mean(rp), ...
            nominal.thd_pct,max(td,[],'omitnan'),mean(td,'omitnan'), ...
            sum(off.thd_valid),min(off.f_sw_Hz),max(off.f_sw_Hz),max(abs(off.speed_error_pct))}; %#ok<AGROW>
    end
end
S=cell2table(rows,'VariableNames',{'param','mode','nominal_ripple_Nm','max_ripple_change_pct', ...
    'mean_ripple_change_pct','nominal_thd_pct','max_thd_change_pct','mean_thd_change_pct', ...
    'valid_off_nominal_thd_count','min_f_sw_Hz','max_f_sw_Hz','max_abs_speed_error_pct'});
S.max_ripple_change_Nm=S.nominal_ripple_Nm.*S.max_ripple_change_pct/100;
S.max_thd_change_percentage_points=S.nominal_thd_pct.*S.max_thd_change_pct/100;
writetable(S,fullfile(out,'sensitivity_summary.csv'));
conv=cell(0,9); fields={'ripple_pp_Nm','thd_pct','f_sw_Hz','speed_mean_rpm'};
for im=1:2
    D=C(C.family=="integration_step"&C.mode==modes(im),:);
    a=D(D.Nsub==20,:); b=D(D.Nsub==40,:);
    for k=1:numel(fields)
        x=a.(fields{k}); y=b.(fields{k}); pct=100*abs(x-y)/max(abs(y),eps);
        conv(end+1,:)={modes(im),string(fields{k}),x,y,pct,5,pct<=5,a.Ts_s/a.Nsub,b.Ts_s/b.Nsub}; %#ok<AGROW>
    end
end
Q=cell2table(conv,'VariableNames',{'mode','metric','N20','N40','relative_change_pct', ...
    'screening_threshold_pct','screening_passed','coarse_step_s','fine_step_s'});
writetable(Q,fullfile(out,'convergence_summary.csv'));
fid=fopen(fullfile(out,'acceptance_summary.txt'),'w','n','UTF-8');
c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Core tests: %d/%d\nMetric tests: %d/%d\nSimulink plant comparisons: %d/%d\n',sum(T.passed),height(T),sum(U.passed),height(U),sum(V.passed),height(V));
fprintf(fid,'Mismatch rows: %d (26 unique configurations; nominal repeated per parameter)\n',height(R));
fprintf(fid,'THD valid: %d/%d; invalid entries retain reason and NaN.\n',sum(R.thd_valid),height(R));
fprintf(fid,'Mean speed gate |error|<2%%: %d/%d\n',sum(abs(R.speed_error_pct)<2),height(R));
fprintf(fid,'Refinement screening <=5%%: %d/%d metrics. This is an engineering screen, not a mathematical convergence proof.\n',sum(Q.screening_passed),height(Q));
fprintf(fid,'No independent experimental replicates; no significance tests or confidence intervals claimed.\n');
for k=1:height(B)
    fprintf(fid,'BASELINE %s speed_rpm=%.6f ripple_Nm=%.6f THD_pct=%.6f fsw_Hz=%.6f recovery_ms=%.6f\n', ...
        B.mode(k),B.speed_mean_rpm(k),B.ripple_pp_Nm(k),B.thd_pct(k),B.f_sw_Hz(k),B.load_recovery_ms(k));
end
disp(B(:,{'mode','speed_mean_rpm','ripple_pp_Nm','thd_pct','f_sw_Hz','load_recovery_ms'}));
disp(S); disp(Q);
end
