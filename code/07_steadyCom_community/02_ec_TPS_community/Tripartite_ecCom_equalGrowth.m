% ============================================================================
% This script saves the equal-growth coupling, defined minimal medium,
% protein-pool bound, and FBA objective directly into the ec community model.

% This allows downstream FBA, FVA, and production-envelope analyses to reproduce
% the equal-growth solution (mu_max = 0.062066/h with
% v_AST = v_PRE = v_POST) without reapplying the constraints all the time.
%
%   1. Equal-flux biomass coupling: v_AST = v_PRE = v_POST.
%   2. Defined minimal medium on the shared community [u] pool.
%   3. Protein-pool bound: AST/PRE/POST_prot_pool_exchange lb = -10 mmol/gDW/h.
%   4. Objective: maximize AST_BIOMASS_reaction.
%
% Input: ecTPS_build.mat   (from Tripartite_ecCom.m)
%
% Output: ecMMU1868-TPS.mat
% ============================================================================

clc; clear;

cd(fileparts(which(mfilename)));

%% Toolbox setup
% Edit toolbox paths
% COBRA Toolbox : https://opencobra.github.io/cobratoolbox/
% Gurobi        : https://www.gurobi.com
COBRA_PATH = 'C:/cobra';
SOLVER_PATH = 'C:/gurobi1202/win64/matlab';

restoredefaultpath; rehash toolboxcache;
addpath(genpath(COBRA_PATH));
addpath(SOLVER_PATH);

initCobraToolbox(false);
changeCobraSolver('gurobi','LP');

%% 1) Input and output paths
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

repoRoot  = fullfile(scriptDir, '..', '..', '..');
modelsDir = fullfile(repoRoot, 'models', '06_community_TPS', 'enzyme_constrained_tps');
srcFile   = fullfile(modelsDir, 'ecTPS_build.mat');
dstFile   = fullfile(modelsDir, 'ecMMU1868-TPS.mat');

fprintf('Loading %s ...\n', srcFile);
S = load(srcFile);
if isfield(S, 'modelCom')
    modelCom = S.modelCom;
else
    f = fieldnames(S); modelCom = S.(f{1});
end
fprintf('  rxns=%d, mets=%d, S nnz=%d\n', ...
    numel(modelCom.rxns), numel(modelCom.mets), nnz(modelCom.S));

% Locate the three biomass reactions
iA = find(strcmp(modelCom.rxns, 'AST_BIOMASS_reaction'),  1);
iP = find(strcmp(modelCom.rxns, 'PRE_BIOMASS_reaction'),  1);
iQ = find(strcmp(modelCom.rxns, 'POST_BIOMASS_reaction'), 1);
assert(~isempty(iA) && ~isempty(iP) && ~isempty(iQ), ...
    'AST_/PRE_/POST_BIOMASS_reaction not found in modelCom.rxns');

%% 2) Equal-flux biomass coupling
% Adds two pseudo-metabolites that force the three biomass reactions to carry
% the same flux: v_AST = v_PRE = v_POST.
nRxn = numel(modelCom.rxns);
row1 = sparse(1, [iA iP], [ 1 -1], 1, nRxn);   %  AST - PRE  = 0
row2 = sparse(1, [iA iQ], [ 1 -1], 1, nRxn);   %  AST - POST = 0
modelCom.S = [modelCom.S; row1; row2];
modelCom.b = [modelCom.b; 0; 0];
if isfield(modelCom, 'csense')
    modelCom.csense = [modelCom.csense(:); 'E'; 'E'];
end
modelCom.mets = [modelCom.mets; {'couple_AST_PRE'; 'couple_AST_POST'}];
modelCom = padMetFields(modelCom, 2);

fprintf('Added equal-flux coupling: mets %d -> %d (+2 pseudo-metabolites)\n', ...
    numel(modelCom.mets) - 2, numel(modelCom.mets));

%% 3) Defined minimal medium (community [u] pool)
% lb/ub pairs from the defined medium used throughout the study.
medium = { ...
  'EX_h2o[u]' -1000  1000;
  'EX_h[u]' -1000  1000;
  'EX_na1[u]' -1000  1000;
  'EX_k[u]' -1000  1000;
  'EX_cl[u]' -1000  1000;
  'EX_ca2[u]' -1000  1000;
  'EX_nh4[u]' -10  1000;
  'EX_pi[u]' -10  1000;
  'EX_so4[u]' -10  1000;
  'EX_o2[u]' -10  1000;
  'EX_co2[u]' 0  1000;
  'EX_glc__D[u]' -10  1000;
  'EX_gln__L[u]' 0  1000;
  'EX_lac__L[u]' 0  1000;
  'EX_hdca[u]' 0  1000;
  'EX_his__L[u]' -0.186207  1000;
  'EX_ile__L[u]' -0.186207  1000;
  'EX_leu__L[u]' -0.186207  1000;
  'EX_lys__L[u]' -0.186207  1000;
  'EX_met__L[u]' -0.155172  1000;
  'EX_phe__L[u]' -0.279310  1000;
  'EX_thr__L[u]' -0.124138  1000;
  'EX_trp__L[u]' -0.034483  1000;
  'EX_val__L[u]' -0.155172  1000 };

% Close every community exchange first
exMask = ~cellfun('isempty', regexp(modelCom.rxns, '^EX_.*\[u\]$'));
exIdx = find(exMask);
modelCom.lb(exIdx) = 0;
modelCom.ub(exIdx) = 1000;
fprintf('Closed %d community EX_*[u] reactions.\n', numel(exIdx));

% Re-open the medium components
missing = {};
for i = 1:size(medium,1)
    rid = medium{i,1};
    j   = find(strcmp(modelCom.rxns, rid), 1);
    if isempty(j)
        missing{end+1} = rid; %#ok<AGROW>
    else
        modelCom.lb(j) = medium{i,2};
        modelCom.ub(j) = medium{i,3};
    end
end
if ~isempty(missing)
    fprintf('medium IDs not in model: %s\n', strjoin(missing, ', '));
end
fprintf('Opened %d/%d defined-medium exchanges.\n', ...
    size(medium,1) - numel(missing), size(medium,1));

%% 4) prot_pool floor (ec only): up to 10 mmol gDW^-1 h^-1 of usable enzyme
% Set prot_pool_exchange lower bound to -10 mmol/gDW/h per cell.
% Sensitivity sweep (-1000 to -0.001) showed growth is constant at 0.1048/h
% from -1000 to -20, drops to 0.0621/h at -10, and becomes infeasible above
% -0.02. The -10 bound gives the physiologically constrained growth rate used here.
for tag = {'AST', 'PRE', 'POST'}
    rid = [tag{1} '_prot_pool_exchange'];
    j = find(strcmp(modelCom.rxns, rid), 1);
    if ~isempty(j)
        modelCom.lb(j) = -10;
    else
        fprintf('WARNING: prot_pool exchange not found: %s\n', rid);
    end
end
fprintf('Set AST/PRE/POST_prot_pool_exchange lb = -10.\n');

%% 5) FBA objective
modelCom.c = zeros(numel(modelCom.rxns), 1);
modelCom.c(iA) = 1;
modelCom.osenseStr = 'max';
fprintf('Objective: c[AST_BIOMASS_reaction] = 1 (osense = max).\n');

%% 6) Verify with optimizeCbModel
fprintf('\nVerifying with optimizeCbModel ...\n');
%% 7) Save BEFORE verify - preserves all annotation fields intact
% optimizeCbModel can strip GECKO annotation fields from the struct,
% so we save first and verify after.
if exist(dstFile, 'file')
    fileattrib(dstFile, '+w', 'a');
    delete(dstFile);
end
save(dstFile, 'modelCom', '-v7.3');
fprintf('\nSaved equal-growth community model: %s\n', dstFile);

%% 8) Verify (after save so field loss does not affect the saved model)
sol = optimizeCbModel(modelCom, 'max');
fprintf(' status = %d\n', sol.stat);
fprintf(' mu_max = %.6f h^-1  \n', sol.f);
fprintf(' AST_BM = %.6f\n', sol.x(iA));
fprintf(' PRE_BM = %.6f\n', sol.x(iP));
fprintf(' POST_BM = %.6f\n', sol.x(iQ));


%% Local Function

function model = padMetFields(model, k)
% Keep metabolite-level fields aligned with model.mets.
    if isfield(model, 'metNames')
        model.metNames  = [model.metNames;  repmat({''}, k, 1)];
    end
    if isfield(model, 'metFormulas')
        model.metFormulas = [model.metFormulas; repmat({''}, k, 1)];
    end
    if isfield(model, 'metCharges')
        model.metCharges = [model.metCharges; zeros(k, 1)];
    end
    if isfield(model, 'metComps')
        model.metComps = [model.metComps; zeros(k, 1)];
    end
end
