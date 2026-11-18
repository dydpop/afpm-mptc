function verify_modeling(repoRoot, outputDir)
% Record standalone port response and repository FOC no-load/load evidence.
here = fileparts(mfilename('fullpath'));
if nargin < 2, outputDir = here; end
if ~isfolder(outputDir), mkdir(outputDir); end
addpath(fullfile(repoRoot,'src','matlab'));
P = afpm_defaults();
for k = 1:2
    Q = P;
    if k == 1, Q.T_load = 0; end
    out = afpm_simulate('foc',[],Q);
    t = out.t;
    rpm = out.wm*60/(2*pi);
    iabc = [out.ia,(-out.ia+sqrt(3)*out.ib)/2,(-out.ia-sqrt(3)*out.ib)/2];
    torque = out.Te;
    steady = t>=1.6 & t<=2;
    time = t(steady);
    angle = out.theta_e(steady);
    coeff = [ones(size(time)),cos(angle),sin(angle)]\iabc(steady,:);
    amplitude = hypot(coeff(2,:),coeff(3,:));
    phase = atan2d(-coeff(3,:),coeff(2,:));
    speedError = 100*(mean(rpm(steady))-P.w_ref)/P.w_ref;
    fprintf('LOAD=%g RPM=%.8f ERROR_PCT=%.8f RMS_A=%s FUND_A=%s PHASE_DEG=%s\n', ...
        Q.T_load,mean(rpm(steady)),speedError,mat2str(rms(iabc(steady,:)),5), ...
        mat2str(amplitude,5),mat2str(phase,5));
    assert(abs(speedError)<2,'Steady speed error exceeds 2%%.');
    label = sprintf('foc_load_%g',Q.T_load);
    writetable(table(t,rpm,iabc(:,1),iabc(:,2),iabc(:,3),torque, ...
        'VariableNames',{'time_s','rpm','ia_A','ib_A','ic_A','torque_Nm'}), ...
        fullfile(outputDir,[label '.csv']));
    fig = figure('Visible','off','Color','w','Position',[0 0 1000 780]);
    layout = tiledlayout(fig,3,1,'TileSpacing','compact');
    title(layout,sprintf('AFPM PWM FOC simulation | load %g N m',Q.T_load));
    ax=nexttile; plot(ax,t,rpm); yline(ax,P.w_ref,'--'); ylabel(ax,'Speed (r/min)'); grid(ax,'on');
    ax=nexttile; plot(ax,t(steady),iabc(steady,:)); xlim(ax,[1.8 1.9]); ylabel(ax,'Phase current (A)'); legend(ax,{'A','B','C'}); grid(ax,'on');
    ax=nexttile; plot(ax,t,torque); ylabel(ax,'Torque (N m)'); xlabel(ax,'Time (s)'); grid(ax,'on');
    exportgraphics(fig,fullfile(outputDir,[label '.png']),'Resolution',150);
    close(fig);
end
modelPath = build_baseline(repoRoot,outputDir);
load_system(modelPath);
cleanup = onCleanup(@() close_system('baseline',0));
dataset = Simulink.SimulationData.Dataset;
dataset = dataset.addElement(timeseries(repmat([1 0 0 1 0 1],2,1),[0;0.001]),'g');
dataset = dataset.addElement(timeseries([0;0],[0;0.001]),'T_load');
dataset = dataset.addElement(timeseries([300;300],[0;0.001]),'w_ref');
input = Simulink.SimulationInput('baseline');
input = input.setExternalInput(dataset);
input = input.setModelParameter('StopTime','0.001','SaveOutput','on','SaveFormat','Dataset');
result = sim(input);
phase = result.yout{1}.Values.Data;
expected = (2*P.Udc/3)/P.Rs*(1-exp(-P.Rs*0.001/P.L))*[1,-0.5,-0.5];
assert(numel(phase(end,:))==3 && max(abs(phase(end,:)-expected))<1e-5, ...
    'Standalone gate-to-current RL check failed.');
assert(max(abs(result.yout{5}.Values.Data),[],'all')<1e-8,'Expected zero electromagnetic torque.');
fprintf('STANDALONE_PASS phase_current_A=%s\n',mat2str(phase(end,:),8));
end
