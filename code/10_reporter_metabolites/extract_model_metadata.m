% This extracts the three model metadata files used by the reporter metabolite
% analysis step. 
%
% Outputs:
%   model_rxn_meta.csv    (rxn_id, rxn_name, subsystem, grRule, lb, ub)
%   model_met_meta.csv    (met_id, met_name, formula)
%   model_S_triplets.csv  (rxn_id, met_id, coefficient)
%
% Run once after generating ecMMU1868-TPS.mat. The output files are
% already provided in this folder.
%
% Input:  ecMMU1868-TPS.mat  (models/06_community_TPS/enzyme_constrained_tps/)
% Output: results/10_reporter_metabolites/ (model_rxn_meta.csv, model_met_meta.csv, model_S_triplets.csv)

%% Toolbox setup
% Edit toolbox paths
%   COBRA Toolbox : https://opencobra.github.io/cobratoolbox/
%   Gurobi        : https://www.gurobi.com
COBRA_PATH  = 'C:/cobra';
SOLVER_PATH = 'C:/gurobi1202/win64/matlab';

restoredefaultpath; rehash toolboxcache;
addpath(genpath(COBRA_PATH));
addpath(SOLVER_PATH);

initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');

%% Paths
here     = pwd;   % run from this folder: cd to code/10_reporter_metabolites 
repoRoot = fullfile(here, '..', '..');

matFile = fullfile(repoRoot, 'models', '06_community_TPS', ...
                   'enzyme_constrained_tps', 'ecMMU1868-TPS.mat');
outDir  = fullfile(repoRoot, 'results', '10_reporter_metabolites');
if ~exist(outDir, 'dir'), mkdir(outDir); end

fprintf('Loading %s ...\n', matFile);
S = load(matFile);
m = S.modelCom;
fprintf('  %d rxns, %d mets, %d genes\n', numel(m.rxns), numel(m.mets), numel(m.genes));

%% Reaction metadata
subSys = strings(numel(m.rxns), 1);
if isfield(m, 'subSystems')
    for j = 1:numel(m.rxns)
        v = m.subSystems{j};
        if iscell(v); v = strjoin(v, ';'); end
        if isstring(v); v = char(v); end
        if isempty(v); v = ''; end
       
        subSys(j) = string(v);
    end
end
Trm = table(m.rxns, m.rxnNames, cellstr(subSys), m.grRules, m.lb, m.ub, ...
    'VariableNames', {'rxn_id','rxn_name','subsystem','grRule','lb','ub'});
writetable(Trm, fullfile(outDir, 'model_rxn_meta.csv'));
fprintf('Wrote model_rxn_meta.csv\n');

%% Metabolite metadata
metForm = repmat({''}, numel(m.mets), 1);
if isfield(m, 'metFormulas'); metForm = m.metFormulas; end
Tmm = table(m.mets, m.metNames, metForm, ...
    'VariableNames', {'met_id','met_name','formula'});
writetable(Tmm, fullfile(outDir, 'model_met_meta.csv'));
fprintf('Wrote model_met_meta.csv\n');

%% S triplets (stoichiometric matrix as rxn_id / met_id / coefficient)
% All nonzero stoichiometric coefficients used to reconstruct the
% metabolite-reaction network for reporter metabolite analysis.
[ii, jj, vv] = find(m.S);
Tst = table(m.rxns(jj), m.mets(ii), vv, ...
    'VariableNames', {'rxn_id','met_id','coef'});
writetable(Tst, fullfile(outDir, 'model_S_triplets.csv'));
fprintf('Wrote model_S_triplets.csv (%d nonzeros)\n', numel(vv));

