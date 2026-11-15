function results = validate_baseline()
% 独立验证 AFPM 等效集中参数 PMSM plant，不验证完整闭环或三维磁场。
% 先给定容差，再运行模型；MATLAB ode45 仅作对照，不向模型注入波形。
modelDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(fileparts(modelDir));
outDir = fullfile(rootDir,'data','research_v2');
if ~isfolder(outDir), mkdir(outDir); end
logPath = fullfile(outDir,'simulink_validation.log');
fid = fopen(logPath,'w','n','UTF-8');
assert(fid >= 0,'无法创建验证日志。');
logCleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'AFPM equivalent lumped-parameter PMSM plant validation only\n');
fprintf(fid,'Not 3D electromagnetic simulation; not full closed-loop validation.\n');
fprintf(fid,'Timestamp=%s\nMATLAB=%s\n',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),version);
fprintf(fid,'SimulinkLicense=%d\n',license('test','Simulink'));
% 在任何仿真前固定判据，各物理量还须分别通过绝对误差判据。
absoluteTolerance = [1e-5,1e-5,1e-5,1e-7,1e-5];
normalizedTolerance = 1e-6;
normalizationFloor = [1,1,1,1,1];
signals = {'ialpha','ibeta','wm','thm','Te'};
units = {'A','A','rad/s','rad','N*m'};
fprintf(fid,'AbsoluteTolerances[A,A,rad/s,rad,N*m]=%s\n',mat2str(absoluteTolerance));
fprintf(fid,'NormalizedTolerance=%.17g; scale=max(reference_peak,1 physical unit)\n',normalizedTolerance);
fprintf(fid,'Simulink: ode45 RelTol=1e-8 AbsTol=1e-10 MaxStep=1e-4 StopTime=0.1\n');
fprintf(fid,'Independent MATLAB: ode45 RelTol=1e-12 AbsTol=1e-13 MaxStep=1e-5\n');
% 将自动生成缓存限制到临时目录，避免修改其他工程文件。
oldConfig = Simulink.fileGenControl('getConfig');
cacheDir = tempname;
mkdir(cacheDir);
Simulink.fileGenControl('set','CacheFolder',cacheDir,'CodeGenFolder',cacheDir);
configCleanup = onCleanup(@() Simulink.fileGenControl('setConfig','config',oldConfig)); %#ok<NASGU>
try
    modelPath = build_baseline_validation();
    model = 'baseline_validation';
    % 从磁盘重新打开，证明模型不依赖构建函数的临时变量。
    close_system(model,0);
    load_system(modelPath);
    modelCleanup = onCleanup(@() close_system(model,0)); %#ok<NASGU>
    integrators = find_system(model,'BlockType','Integrator');
    fromWorkspace = find_system(model,'BlockType','FromWorkspace');
    assert(numel(integrators)==4,'应具有四个独立连续积分状态。');
    assert(isempty(fromWorkspace),'禁止用 From Workspace 输入预计算状态。');
    fprintf(fid,'ReloadedModel=%s\nIntegratorCount=%d; FromWorkspaceCount=%d\n',modelPath,numel(integrators),numel(fromWorkspace));
    % 参考侧独立编码参数和方程，不调用模型 Fcn 表达式。
    P = struct('Rs',0.2,'L',0.0085,'psi',0.175,'p',4,'J',0.003,'F',0.008);
    fprintf(fid,'Rs=0.2 L=0.0085 psi=0.175 p=4 J=0.003 F=0.008; x0=[0 0 0 0]; Tload=0\n');
    rows = cell(0,11);
    cases = {'rotating_voltage_20Hz','aligned_DC_RL_analytic'};
    for testIndex = 1:numel(cases)
        input = Simulink.SimulationInput(model);
        frequency = 20;
        if testIndex == 2
            % DC 电压对齐 alpha 轴，thm=ibeta=wm=0 是不变子空间。
            % 因 Te=0，未额外钳制机械状态，退化为严格可解析的 RL 测试。
            frequency = 0;
            input = input.setVariable('inputFrequency',0,'Workspace',model);
        end
        simOut = sim(input);
        trace = simOut.get('plant_trace');
        t = double(trace.Time(:));
        actual = double(squeeze(trace.Data));
        assert(size(actual,1)==numel(t) && size(actual,2)==5,'输出维度必须为 N×5。');
        assert(numel(t)>100 && all(diff(t)>0),'求解器时间点数量或顺序不正确。');
        assert(abs(t(1))<1e-14 && abs(t(end)-0.1)<1e-12,'仿真没有完整运行 0.1 秒。');
        assert(all(isfinite(actual),'all'),'Simulink 输出含非有限值。');
        options = odeset('RelTol',1e-12,'AbsTol',1e-13,'MaxStep',1e-5);
        referenceSolution = ode45(@(time,x) reference_rhs(time,x,P,frequency),[0 0.1],zeros(4,1),options);
        states = deval(referenceSolution,t).';
        thetaElectrical = P.p*states(:,4);
        torque = 1.5*P.p*P.psi*(-states(:,1).*sin(thetaElectrical)+states(:,2).*cos(thetaElectrical));
        reference = [states torque];
        fprintf(fid,'Test=%s Samples=%d FinalTime=%.17g\n',cases{testIndex},numel(t),t(end));
        fprintf(fid,'FinalStateAndTorque=%s\n',mat2str(actual(end,:),12));
        comparisons = {reference};
        methods = {'independent_ode45'};
        if testIndex == 2
            analytic = zeros(numel(t),5);
            analytic(:,1) = (30/P.Rs)*(1-exp(-P.Rs*t/P.L));
            comparisons{end+1} = analytic; %#ok<AGROW>
            methods{end+1} = 'exact_RL_solution'; %#ok<AGROW>
        end
        for methodIndex = 1:numel(comparisons)
            expected = comparisons{methodIndex};
            errors = max(abs(actual-expected),[],1);
            scale = max(max(abs(expected),[],1),normalizationFloor);
            normalized = errors./scale;
            passed = errors <= absoluteTolerance & normalized <= normalizedTolerance;
            for k = 1:5
                rows(end+1,:) = {cases{testIndex},methods{methodIndex},signals{k},units{k}, ...
                    errors(k),normalized(k),scale(k),absoluteTolerance(k), ...
                    normalizedTolerance,passed(k),numel(t)}; %#ok<AGROW>
                fprintf(fid,'%s %s max_abs=%.12g normalized=%.12g pass=%d\n', ...
                    methods{methodIndex},signals{k},errors(k),normalized(k),passed(k));
            end
        end
    end
    results = cell2table(rows,'VariableNames',{'test','reference','signal','unit', ...
        'max_abs_error','normalized_error','normalization_scale', ...
        'absolute_tolerance','normalized_tolerance','passed','sample_count'});
    csvPath = fullfile(outDir,'simulink_validation.csv');
    writetable(results,csvPath);
    readback = readtable(csvPath);
    assert(height(readback)==15 && all(readback.passed),'保存的验证结果存在失败或缺行。');
    fprintf(fid,'CSV=%s\nOVERALL_PASS=%d\n',csvPath,all(results.passed));
    disp(results);
    fprintf('VALIDATION_PASS=%d\nCSV=%s\nLOG=%s\n',all(results.passed),csvPath,logPath);
    assert(all(results.passed),'独立 plant 验证未通过预设阈值。');
catch err
    fprintf(fid,'FAILED: %s\n',getReport(err,'extended','hyperlinks','off'));
    rethrow(err);
end
end

function dx = reference_rhs(t,x,P,frequency)
% 独立 MATLAB 参考模型：旋转电压、反电动势和机械方程显式展开。
ua = 30*cos(2*pi*frequency*t);
ub = 30*sin(2*pi*frequency*t);
thetaElectrical = P.p*x(4);
omegaElectrical = P.p*x(3);
ia = x(1);
ib = x(2);
te = 1.5*P.p*P.psi*(-ia*sin(thetaElectrical)+ib*cos(thetaElectrical));
dx = [(ua-P.Rs*ia+P.psi*omegaElectrical*sin(thetaElectrical))/P.L; ...
      (ub-P.Rs*ib-P.psi*omegaElectrical*cos(thetaElectrical))/P.L; ...
      (te-P.F*x(3))/P.J; x(3)];
end
