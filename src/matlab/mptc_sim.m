function varargout = mptc_sim(task)
% Compatibility entry for the corrected, tested MATLAB research workflow.
% The original ideal-voltage FOC and inconsistent delay implementation are
% preserved in Git history, not retained as a second active experiment.
% demo: regression tests + matched nominal baseline; no incomplete-batch plots.
% sweep: full validated batch including mismatch, robustness and plots.
if nargin<1, task='demo'; end
switch lower(task)
    case 'demo'
        run_research('tests'); R=run_research('baseline');
    case 'sweep'
        R=run_research('all');
    otherwise
        R=run_research(task);
end
if nargout>0, varargout{1}=R; end
end
