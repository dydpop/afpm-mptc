function modelPath = build_baseline_validation()
% AFPM 等效集中参数 PMSM 验证对象；不是三维磁场仿真。
% 连续四状态完全由 Simulink Integrator 求解，禁止预计算状态波形输入。
% 状态顺序：[ialpha, ibeta, wm, thm]，角度为机械角度。
modelDir = fileparts(mfilename('fullpath'));
model = 'baseline_validation';
modelPath = fullfile(modelDir, [model '.slx']);
assert(license('test','Simulink') == 1, '需要可用的 Simulink 许可。');
if bdIsLoaded(model)
    assert(strcmp(get_param(model,'Dirty'),'off'), '模型存在未保存修改，请先保存。');
    close_system(model,0);
end
load_system('simulink');
new_system(model);
set_param(model,'SolverType','Variable-step','Solver','ode45', ...
    'RelTol','1e-8','AbsTol','1e-10','MaxStep','1e-4', ...
    'StartTime','0','StopTime','0.1','ReturnWorkspaceOutputs','on', ...
    'SaveTime','on','TimeSaveName','tout','SaveOutput','off', ...
    'SignalLogging','off');
set_param(model,'Description', ...
    'AFPM 等效集中参数、非凸极 PMSM plant 独立验证；不代表三维磁场或完整闭环控制。');
w = get_param(model,'ModelWorkspace');
params = struct('Rs',0.2,'L',0.0085,'psi',0.175,'p',4, ...
    'J',0.003,'F',0.008,'voltageAmplitude',30,'inputFrequency',20,'loadTorque',0);
names = fieldnames(params);
for k = 1:numel(names)
    assignin(w,names{k},params.(names{k}));
end
add_block('simulink/Sources/Sine Wave',[model '/ua_V'], ...
    'Amplitude','voltageAmplitude','Frequency','2*pi*inputFrequency', ...
    'Phase','pi/2','SampleTime','0','Position',[35 55 95 90]);
add_block('simulink/Sources/Sine Wave',[model '/ub_V'], ...
    'Amplitude','voltageAmplitude','Frequency','2*pi*inputFrequency', ...
    'Phase','0','SampleTime','0','Position',[35 130 95 165]);
add_block('simulink/Sources/Constant',[model '/Tload_Nm'], ...
    'Value','loadTorque','Position',[35 205 95 240]);
plant = [model '/AFPM_equivalent_PMSM'];
add_block('simulink/Ports & Subsystems/Subsystem',plant, ...
    'Position',[195 45 420 265]);
Simulink.SubSystem.deleteContents(plant);
inNames = {'ua','ub','Tload'};
inUnits = {'V','V','N*m'};
for k = 1:3
    add_block('simulink/Sources/In1',[plant '/' inNames{k}], ...
        'Port',num2str(k),'Unit',inUnits{k}, ...
        'Position',[25 405+55*k 55 425+55*k]);
end
stateNames = {'ialpha','ibeta','wm','thm'};
stateUnits = {'A','A','rad/s','rad'};
add_block('simulink/Signal Routing/Mux',[plant '/state_and_input'], ...
    'Inputs','7','Position',[185 50 190 400]);
expressions = { ...
    '(u(5)-Rs*u(1)+psi*p*u(3)*sin(p*u(4)))/L', ...
    '(u(6)-Rs*u(2)-psi*p*u(3)*cos(p*u(4)))/L', ...
    '(1.5*p*psi*(-u(1)*sin(p*u(4))+u(2)*cos(p*u(4)))-u(7)-F*u(3))/J', ...
    'u(3)'};
for k = 1:4
    y = 35 + 100*(k-1);
    add_block('simulink/User-Defined Functions/Fcn',[plant '/d_' stateNames{k}], ...
        'Expr',expressions{k},'Position',[285 y 465 y+40]);
    add_block('simulink/Continuous/Integrator',[plant '/integrate_' stateNames{k}], ...
        'InitialCondition','0','Position',[520 y 555 y+40]);
    add_block('simulink/Sinks/Out1',[plant '/' stateNames{k}], ...
        'Port',num2str(k),'Unit',stateUnits{k}, ...
        'Position',[650 y+10 680 y+30]);
    add_line(plant,'state_and_input/1',['d_' stateNames{k} '/1'],'autorouting','on');
    add_line(plant,['d_' stateNames{k} '/1'],['integrate_' stateNames{k} '/1'],'autorouting','on');
    line = add_line(plant,['integrate_' stateNames{k} '/1'],[stateNames{k} '/1'],'autorouting','on');
    set_param(line,'Name',stateNames{k});
    add_line(plant,['integrate_' stateNames{k} '/1'],['state_and_input/' num2str(k)],'autorouting','on');
end
for k = 1:3
    line = add_line(plant,[inNames{k} '/1'],['state_and_input/' num2str(k+4)],'autorouting','on');
    set_param(line,'Name',inNames{k});
end
add_block('simulink/User-Defined Functions/Fcn',[plant '/electromagnetic_torque'], ...
    'Expr','1.5*p*psi*(-u(1)*sin(p*u(4))+u(2)*cos(p*u(4)))', ...
    'Position',[285 450 465 495]);
add_block('simulink/Sinks/Out1',[plant '/Te'], ...
    'Port','5','Unit','N*m','Position',[650 460 680 480]);
add_line(plant,'state_and_input/1','electromagnetic_torque/1','autorouting','on');
line = add_line(plant,'electromagnetic_torque/1','Te/1','autorouting','on');
set_param(line,'Name','Te');
inputBlocks = {'ua_V','ub_V','Tload_Nm'};
for k = 1:3
    line = add_line(model,[inputBlocks{k} '/1'],['AFPM_equivalent_PMSM/' num2str(k)],'autorouting','on');
    set_param(line,'Name',inNames{k});
end
add_block('simulink/Signal Routing/Mux',[model '/trace_mux'], ...
    'Inputs','5','Position',[520 55 525 265]);
outputNames = [stateNames {'Te'}];
for k = 1:5
    line = add_line(model,['AFPM_equivalent_PMSM/' num2str(k)], ...
        ['trace_mux/' num2str(k)],'autorouting','on');
    set_param(line,'Name',outputNames{k});
end
add_block('simulink/Sinks/To Workspace',[model '/solver_output'], ...
    'VariableName','plant_trace','SaveFormat','Timeseries', ...
    'MaxDataPoints','inf','Decimation','1', ...
    'Position',[610 125 735 170]);
add_line(model,'trace_mux/1','solver_output/1','autorouting','on');
annotation = Simulink.Annotation(model, ...
    sprintf('AFPM 等效集中参数 PMSM（非三维磁场）\n四状态连续积分；30 V / 20 Hz 旋转电压；零初值、零负载\n输出：ialpha [A], ibeta [A], wm [rad/s], thm [rad], Te [N*m]'));
annotation.Position = [30 300 750 365];
set_param(model,'SimulationCommand','update');
save_system(model,modelPath);
fprintf('MODEL_SAVED=%s\n',modelPath);
end
