function modelPath = build_baseline(repoRoot, outputDir)
% Build the standalone AFPM plant using the repository's continuous solver.
here = fileparts(mfilename('fullpath'));
if nargin < 2, outputDir = here; end
if ~isfolder(outputDir), mkdir(outputDir); end
addpath(fullfile(repoRoot,'src','matlab'));
P = afpm_defaults();
model = 'baseline';
modelPath = fullfile(outputDir,[model '.slx']);
assert(~bdIsLoaded(model),'Close the existing baseline model before building.');
source = 'baseline_validation';
wasLoaded = bdIsLoaded(source);
load_system(fullfile(repoRoot,'src','model',[source '.slx']));
if ~wasLoaded, sourceCleanup = onCleanup(@() close_system(source,0)); end %#ok<NASGU>
load_system('simulink');
new_system(model);
cleanup = onCleanup(@() close_system(model,0));
w = get_param(model,'ModelWorkspace');
for name = {'Rs','L','psi','p','J','F','Udc'}
    assignin(w,name{1},P.(name{1}));
end
set_param(model,'Solver','ode45','RelTol','1e-8','AbsTol','1e-10', ...
    'MaxStep','1e-5','StopTime','0.1','ReturnWorkspaceOutputs','on');
plant = 'AFPM_equivalent_PMSM';
add_block([source '/' plant],[model '/' plant],'Position',[365 53 505 190]);
inNames = {'g','T_load','w_ref'};
widths = [6 1 1];
inUnits = {'1','N*m','rpm'};
for k=1:3
    add_block('simulink/Sources/In1',[model '/' inNames{k}], ...
        'Port',num2str(k),'PortDimensions',num2str(widths(k)), ...
        'Unit',inUnits{k},'Position',[25 30+75*k 55 44+75*k]);
end
add_block('simulink/Sinks/Terminator',[model '/Test_reference_only'], ...
    'Position',[110 254 130 270]);
add_line(model,'w_ref/1','Test_reference_only/1','autorouting','on');
voltage = {'Udc*(2*(u(1)-u(2))-(u(3)-u(4))-(u(5)-u(6)))/6', ...
    'Udc*((u(3)-u(4))-(u(5)-u(6)))/(2*sqrt(3))'};
for k=1:2
    name = sprintf('Voltage_%d',k);
    add_block('simulink/User-Defined Functions/Fcn',[model '/' name], ...
        'Expr',voltage{k},'Position',[170 18+45*k 292 48+45*k]);
    add_line(model,'g/1',[name '/1'],'autorouting','on');
    add_line(model,[name '/1'],sprintf('%s/%d',plant,k),'autorouting','on');
end
add_line(model,'T_load/1',[plant '/3'],'autorouting','on');
add_block('simulink/Signal Routing/Mux',[model '/State'],'Inputs','4', ...
    'Position',[550 40 555 188]);
for k=1:4
    add_line(model,sprintf('%s/%d',plant,k),sprintf('State/%d',k),'autorouting','on');
end
outNames = {'i_abc','idq','w_r','theta_e','T_e','psi_ab','Vdc'};
units = {'A','A','rad/s','rad','N*m','Wb','V'};
expressions = {{'u(1)','(-u(1)+sqrt(3)*u(2))/2','(-u(1)-sqrt(3)*u(2))/2'}, ...
    {'u(1)*cos(p*u(4))+u(2)*sin(p*u(4))','-u(1)*sin(p*u(4))+u(2)*cos(p*u(4))'}, ...
    {}, {'p*u(4)'}, {}, {'L*u(1)+psi*cos(p*u(4))','L*u(2)+psi*sin(p*u(4))'}, {}};
for k=1:7
    y = 35+115*(k-1);
    add_block('simulink/Sinks/Out1',[model '/' outNames{k}], ...
        'Port',num2str(k),'Unit',units{k},'Position',[1020 y 1050 y+16]);
    parts = expressions{k};
    switch outNames{k}
        case 'w_r'
            output = [plant '/3'];
        case 'T_e'
            output = [plant '/5'];
        case 'Vdc'
            add_block('simulink/Sources/Constant',[model '/DC_bus'], ...
                'Value','Udc','Position',[860 y 910 y+19]);
            output = 'DC_bus/1';
        otherwise
            if numel(parts)>1
                mux = [outNames{k} '_vector'];
                add_block('simulink/Signal Routing/Mux',[model '/' mux], ...
                    'Inputs',num2str(numel(parts)),'Position',[945 y 950 y+85]);
                output = [mux '/1'];
            end
            for j=1:numel(parts)
                name = sprintf('%s_component_%d',outNames{k},j);
                add_block('simulink/User-Defined Functions/Fcn',[model '/' name], ...
                    'Expr',parts{j},'Position',[670 y+31*(j-1) 880 y+22+31*(j-1)]);
                add_line(model,'State/1',[name '/1'],'autorouting','on');
                if numel(parts)>1
                    add_line(model,[name '/1'],sprintf('%s/%d',mux,j),'autorouting','on');
                else
                    output = [name '/1'];
                end
            end
    end
    add_line(model,output,[outNames{k} '/1'],'autorouting','on');
end
set_param(model,'Description',['Lumped AFPM equivalent with ideal complementary gates ' ...
    '[A+ A- B+ B- C+ C-]. w_ref is a test interface only.']);
set_param(model,'SimulationCommand','update');
save_system(model,modelPath);
end
