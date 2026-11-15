function R = run_research(task)
% Reproducible research stage. All motor simulations execute in MATLAB.
% Tests, baseline, sweep, robustness, figures, all. No claim is predetermined.
if nargin < 1, task = 'all'; end
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'src','matlab'));
out = fullfile(root,'data','research_v2');
figdir = fullfile(root,'figures','research_v2');
if ~isfolder(out), mkdir(out); end
if ~isfolder(figdir), mkdir(figdir); end
P = afpm_defaults(); P.metric_window_s = 0.4;
fingerprint = source_fingerprint(root);
manifestFile = fullfile(out,'manifest.mat');
if any(strcmp(task,{'all','baseline'}))
    manifest = struct('run_id',['run_' char(datetime('now','Format','yyyyMMdd_HHmmss'))], ...
        'source_fingerprint',fingerprint,'status','running','expected_sweep_rows',30, ...
        'expected_robustness_rows',17,'completed_utc_note','local time recorded');
    save(manifestFile,'manifest');
    writetable(struct2table(manifest),fullfile(out,'manifest.csv'));
elseif any(strcmp(task,{'sweep','robustness','figures'}))
    assert(isfile(manifestFile),'Run all first to establish an auditable batch.');
    v=load(manifestFile,'manifest');
    assert(strcmp(v.manifest.source_fingerprint,fingerprint),'Source changed since this batch. Rerun all.');
end
P.source_fingerprint=fingerprint;
if any(strcmp(task,{'all','baseline'})), P.run_id=manifest.run_id;
else, P.run_id=['run_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]; end
R = table();
if any(strcmp(task,{'tests','all'}))
    A = afpm_core_tests(); B = afpm_metric_tests();
    if isstruct(A), A = struct2table(A); end
    if ismember('pass',A.Properties.VariableNames), A.Properties.VariableNames{'pass'} = 'passed'; end
    if isstruct(B), B = struct2table(B); end
    writetable(A,fullfile(out,'core_tests.csv'));
    writetable(B,fullfile(out,'metric_tests.csv'));
    disp(A); disp(B);
    assert(all(A.passed) && all(B.passed),'Research:TestsFailed','Regression tests failed.');
end
if any(strcmp(task,{'baseline','all'}))
    diary(fullfile(out,'baseline_log.txt')); clean = onCleanup(@() diary('off'));
    fprintf('MATLAB %s | %s\n',version, char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
    o = afpm_simulate('mptc',struct(),P);
    m = afpm_metrics(o,P);
    % Match measured nominal switching frequency, not a guessed fixed value.
    assert(isfinite(m.f_sw_Hz) && m.f_sw_Hz>0);
    P.f_sw_foc = m.f_sw_Hz;
    [R,~] = add_case(R,'baseline_mptc','baseline','mptc','none',1,P,o,out);
    of = afpm_simulate('foc',struct(),P);
    [R,~] = add_case(R,'baseline_foc','baseline','foc','none',1,P,of,out);
    save(fullfile(out,'protocol.mat'),'P','-v7');
    writetable(R,fullfile(out,'baseline_metrics.csv'));
    protocol = struct2table(P,'AsArray',true);
    writetable(protocol,fullfile(out,'protocol.csv'));
    disp(R(:,{'case_id','mode','speed_mean_rpm','ripple_pp_Nm','thd_pct','f_sw_Hz'}));
    fprintf('Nominal carrier calibrated once to %.6f Hz; held fixed in all mismatch cases.\n',P.f_sw_foc);
    clear clean;
end
if any(strcmp(task,{'sweep','robustness','all'}))
    assert(isfile(fullfile(out,'protocol.mat')),'Run baseline first.');
    s = load(fullfile(out,'protocol.mat'),'P'); P = s.P;
end
if any(strcmp(task,{'sweep','all'}))
    R = table(); modes = {'mptc','foc'}; pars = {'Rs','L','psi'}; factors = [0.8 0.9 1 1.1 1.2];
    for im = 1:2
        for ip = 1:3
            for ix = 1:5
                mm = struct(pars{ip},factors(ix));
                id = sprintf('%s_%s_%03d',modes{im},pars{ip},round(100*factors(ix)));
                fprintf('SWEEP %s\n',id);
                o = afpm_simulate(modes{im},mm,P);
                [R,~] = add_case(R,id,'mismatch',modes{im},pars{ip},factors(ix),P,o,out);
                writetable(R,fullfile(out,'sweep_metrics.csv'));
            end
        end
    end
    assert(height(R)==30); disp(R);
end
if any(strcmp(task,{'robustness','all'}))
    R = table();
    for mode = {'mptc','foc'}
        for nsub = [10 20 40]
            Q = P; Q.Nsub = nsub;
            id = sprintf('step_%s_n%02d',mode{1},nsub);
            fprintf('CONVERGENCE %s\n',id);
            o = afpm_simulate(mode{1},struct(),Q);
            [R,~] = add_case(R,id,'integration_step',mode{1},'Nsub',nsub,Q,o,out);
            writetable(R,fullfile(out,'robustness_metrics.csv'));
        end
    end
    for ts = [50e-6 100e-6]
        for enabled = [false true]
            Q = P; Q.Ts = ts; Q.Nsub = round(ts/(P.Ts/P.Nsub)); Q.delay_comp = enabled;
            id = sprintf('delay_%dus_%d',round(ts*1e6),enabled);
            o = afpm_simulate('mptc',struct(),Q);
            [R,~] = add_case(R,id,'delay_ablation','mptc','delay_comp',double(enabled),Q,o,out);
            writetable(R,fullfile(out,'robustness_metrics.csv'));
        end
    end
    for weight = [0 P.lam_sw 0.05]
        Q = P; Q.lam_sw = weight;
        id = sprintf('switch_weight_%03d',round(1000*weight));
        o = afpm_simulate('mptc',struct(),Q);
        [R,~] = add_case(R,id,'switch_weight','mptc','lam_sw',weight,Q,o,out);
        writetable(R,fullfile(out,'robustness_metrics.csv'));
    end
    for weight = [0.05 0.5 1 2]
        Q = P; Q.lam_flux = weight;
        id = sprintf('flux_weight_%03d',round(100*weight));
        o = afpm_simulate('mptc',struct(),Q);
        [R,~] = add_case(R,id,'flux_weight','mptc','lam_flux',weight,Q,o,out);
        writetable(R,fullfile(out,'robustness_metrics.csv'));
    end
    disp(R(:,{'case_id','mode','ripple_pp_Nm','thd_pct','f_sw_Hz'}));
end
if any(strcmp(task,{'figures','all'}))
    A=readtable(fullfile(out,'sweep_metrics.csv'),'TextType','string');
    B=readtable(fullfile(out,'robustness_metrics.csv'),'TextType','string');
    assert(height(A)==30 && height(B)==17,'Incomplete batch: do not publish figures.');
    v=load(manifestFile,'manifest');
    assert(all(A.source_fingerprint==string(fingerprint)) && all(B.source_fingerprint==string(fingerprint)),'Mixed source versions.');
    assert(all(A.run_id==string(v.manifest.run_id)) && all(B.run_id==string(v.manifest.run_id)),'Mixed batches.');
    assert(strcmp(source_fingerprint(root),fingerprint),'Sources changed during execution.');
    make_figures(out,figdir);
end
if any(strcmp(task,{'all','figures'}))
    v=load(manifestFile,'manifest'); manifest=v.manifest;
    manifest.status='complete';
    manifest.completed_local=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    save(manifestFile,'manifest');
    writetable(struct2table(manifest),fullfile(out,'manifest.csv'));
end
end

function fingerprint=source_fingerprint(root)
files={'afpm_defaults.m','afpm_simulate.m','afpm_metrics.m','afpm_core_tests.m','afpm_metric_tests.m','run_research.m'};
md=java.security.MessageDigest.getInstance('SHA-256');
for k=1:numel(files)
    bytes=unicode2native([files{k} char(10) fileread(fullfile(root,'src','matlab',files{k}))],'UTF-8');
    md.update(typecast(uint8(bytes),'int8'));
end
fingerprint=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end

function [R,m] = add_case(R,id,family,mode,param,factor,P,o,out)
m = afpm_metrics(o,P);
q = o.t>=P.t_end-P.metric_window_s;
m.id_mean_A = mean(o.id(q)); m.id_control_peak_abs_A = max(abs(o.id(q)));
m.phaseA_control_peak_A = max(abs(o.ia(q)));
m.peak_sample_source = 'control-rate phaseA/id, not all-phase or analog peak';
m.foc_Kp_if_used = o.controller_params.L*400;
m.foc_Ki_if_used = o.controller_params.Rs*400;
m.parameter_protocol = 'estimated parameters configure both FOC gains and feedforward; MPTC weights fixed';
m.run_id = P.run_id; m.source_fingerprint = P.source_fingerprint;
m.flux_mean_Wb = mean(hypot(P.psi+P.L*o.id(q),P.L*o.iq(q)));
m.controller_Rs = o.controller_params.Rs; m.controller_L = o.controller_params.L;
m.controller_psi = o.controller_params.psi;
base = table(string(id),string(family),string(mode),string(param),factor,P.Ts,P.Nsub,P.f_sw_foc,...
    'VariableNames',{'case_id','family','mode','param','factor','Ts_s','Nsub','carrier_Hz'});
row = [base struct2table(m,'AsArray',true)];
R = [R; row];
% Complete trajectories stay local; CSV summaries and selected waveforms are versioned.
save(fullfile(out,[id '.mat']),'P','o','m','-v7');
if strcmp(family,'baseline')
    W = table(o.t(:),o.wm(:)*60/(2*pi),o.Te(:),o.Tref(:),o.ia(:),o.ib(:),o.edges(:),...
        'VariableNames',{'time_s','speed_rpm','torque_Nm','torque_ref_Nm','i_alpha_A','i_beta_A','bridge_edges'});
    writetable(W,fullfile(out,[id '_waveform.csv']));
end
fprintf('%s | speed %.3f rpm | ripple %.5f Nm | THD %.5f %% | fsw %.3f Hz\n',...
    id,m.speed_mean_rpm,m.ripple_pp_Nm,m.thd_pct,m.f_sw_Hz);
end

function make_figures(out,figdir)
M = load(fullfile(out,'baseline_mptc.mat')); F = load(fullfile(out,'baseline_foc.mat'));
colors = [0 0.4470 0.7410;0.8500 0.3250 0.0980];
f = figure('Visible','off','Color','w','Position',[80 80 1100 800]);
tiledlayout(3,1,'TileSpacing','compact');
nexttile; plot(M.o.t,M.o.wm*60/(2*pi),'Color',colors(1,:)); hold on;
plot(F.o.t,F.o.wm*60/(2*pi),'Color',colors(2,:)); yline(M.P.w_ref,'k--');
xline(M.P.t_step,':'); ylabel('Speed (rpm)'); grid on; legend('MPTC','PWM FOC','Reference','Location','best');
title('Same plant, torque outer loop, bus and one-sample actuation delay');
nexttile; plot(M.o.t,M.o.Te,'Color',colors(1,:)); hold on; plot(F.o.t,F.o.Te,'Color',colors(2,:));
xlim([M.P.t_end-0.02 M.P.t_end]); ylabel('Torque (N m)'); grid on; title('Steady-state switching ripple (control-rate display)');
nexttile; plot(M.o.t,M.o.ia,'Color',colors(1,:)); hold on; plot(F.o.t,F.o.ia,'Color',colors(2,:));
xlim([M.P.t_end-0.15 M.P.t_end]); ylabel('Phase-a current (A)'); xlabel('Time (s)'); grid on;
exportgraphics(f,fullfile(figdir,'baseline_comparison.png'),'Resolution',300); close(f);
if isfile(fullfile(out,'sweep_metrics.csv'))
    R = readtable(fullfile(out,'sweep_metrics.csv'),'TextType','string');
    f = figure('Visible','off','Color','w','Position',[80 80 1150 1100]); tiledlayout(4,3,'TileSpacing','compact');
    pars = ["Rs","L","psi"]; modes = ["mptc","foc"];
    yfields={'ripple_pp_Nm','thd_pct','load_recovery_ms','f_sw_Hz'};
    ylabels={'Torque ripple p-p (N m)','THD, harmonics 2-40 (%)','Load speed recovery (ms)','Device switching frequency (Hz)'};
    for metric = 1:4
        for ip = 1:3
            nexttile; hold on;
            for im = 1:2
                D = R(R.param==pars(ip)&R.mode==modes(im),:);
                y = D.(yfields{metric});
                plot((D.factor-1)*100,y,'-o','Color',colors(im,:),'LineWidth',1.3);
            end
            grid on; xlabel('Controller parameter error (%)'); title(pars(ip),'Interpreter','none');
            ylabel(ylabels{metric});
            if ip==1, legend('MPTC','PWM FOC','Location','best'); end
        end
    end
    sgtitle('Fixed physical plant; nominally calibrated PWM carrier held fixed');
    exportgraphics(f,fullfile(figdir,'mismatch_sensitivity.png'),'Resolution',300); close(f);
end
if isfile(fullfile(out,'robustness_metrics.csv'))
    R = readtable(fullfile(out,'robustness_metrics.csv'),'TextType','string');
    f = figure('Visible','off','Color','w','Position',[80 80 1000 400]); tiledlayout(1,2);
    nexttile; hold on; modes = ["mptc","foc"];
    for im=1:2
        D=R(R.family=="integration_step"&R.mode==modes(im),:);
        plot(D.Ts_s./D.Nsub*1e6,D.ripple_pp_Nm,'-o','Color',colors(im,:));
    end
    grid on; xlabel('Plant integration step (us)'); ylabel('Torque ripple, p-p (N m)'); legend('MPTC','PWM FOC');
    title('Numerical refinement: smaller step must be checked');
    nexttile; D=R(R.family=="switch_weight",:);
    plot(D.f_sw_Hz,D.ripple_pp_Nm,'-o'); grid on; xlabel('Measured device switching frequency (Hz)');
    ylabel('Torque ripple, p-p (N m)'); title('Switching-penalty ablation (not an optimized Pareto front)');
    exportgraphics(f,fullfile(figdir,'robustness_checks.png'),'Resolution',300); close(f);
end
end
