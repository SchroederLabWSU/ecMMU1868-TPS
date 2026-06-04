% ============================================================================
% GECKO 3.0 ecModel reconstruction protocol for the astrocyte+PAP(iMMU1868a-pap) model.
%
% Stages:
%   0  Installation & adapter setup
%   1  Conventional GEM expansion to ecModel structure -> ecModel (UniProt + ComplexPortal)
%   2  Kcat integration (BRENDA + DLKcat) into ecModel
%   3  Sensitivity tuning to recover max growth
%   4  Proteomics integration and flexibilization to gR_exp
%   5  Validation, enzyme usage, ecModel vs conventional flux comparison
%
% The YMLs model outputs (stage1-4) and the final YML ecModel are written to
% ./models/. The tuned kcats and run log are written to ./output/.
% A final .mat model (for downstream SteadyCom) is also written to ./output/.
% ============================================================================

%% Stage 0 - Installation
% startGECKOproject() is intentionally omitted. It is an interactive
% evironment that generates a new adapter and folder structure, which 
% would overwrite the existing configured adapter and model files already here.
checkInstallation
GECKOInstaller.install
setRavenSolver('gurobi')

% Locate this protocol's adapter
adapterLocation = fullfile(fileparts(mfilename('fullpath')), ...
                           'ecMMU1868aPapAdapter.m');
ModelAdapterManager.setDefault(adapterLocation);
params = ModelAdapterManager.getDefault().getParameters();

%% Stage 1 - Build ecModel structure using the conventional GEM
model = loadConventionalGEM();

assert(any(strcmp(model.rxns, params.c_source)), ...
    'Carbon-source reaction %s not found in iMMU1868a-pap.', params.c_source);
assert(any(strcmp(model.rxns, params.bioRxn)), ...
    'Biomass reaction %s not found in iMMU1868a-pap.', params.bioRxn);

[ecModel, noUniprot] = makeEcModel(model, false);
fprintf('Genes without UniProt mapping: %d / %d\n', ...
    numel(noUniprot), numel(model.genes));

[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);
saveEcModel(ecModel, 'ecMMU1868a-pap_stage1.yml');

coverage = (numel(model.genes) - numel(noUniprot)) / numel(model.genes) * 100;
fprintf('UniProt gene coverage: %.2f %%\n', coverage);

%% Stage 2 - Kcat integration
ecModel = getECfromGEM(ecModel);
ecModel = getECfromDatabase(ecModel);

kcatList_fuzzy = fuzzyKcatMatching(ecModel);
[ecModel, noSmiles] = findMetSmiles(ecModel);

% writeDLKcatInput(ecModel,[],[],[],[],true);  % uncomment to regenerate
% runDLKcat();                                  % requires Docker

kcatList_DLKcat = readDLKcatOutput(ecModel); % only if you already have DLKcat.tsv
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy);

ecModel = selectKcatValue(ecModel, kcatList_merged);
ecModel = getKcatAcrossIsozymes(ecModel);
[ecModel, rxnsMissingGPR, standardMW, standardKcat] = getStandardKcat(ecModel);
ecModel = applyKcatConstraints(ecModel);

f = params.f;   % picks value from adapter file, recomputed in Stage 4 from proteomics
ecModel = setProtPoolSize(ecModel, params.Ptot, f, params.sigma);
saveEcModel(ecModel, 'ecMMU1868a-pap_stage2.yml');

%% Sanity check - kcat-constrained model growth check at unlimited glucose
ecModel = setParam(ecModel, 'lb',  params.c_source, -1000);
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);
sol = solveLP(ecModel, 1);
fprintf('Stage 2 ecModel growth (glc unlimited): %.4f 1/h\n', sol.f);

% Reference with conventional GEM growth at the same glucose setting
stoic_model = loadConventionalGEM();
stoic_model = setParam(stoic_model, 'lb',  params.c_source, -1000);
stoic_model = setParam(stoic_model, 'obj', params.bioRxn,    1);
soln = solveLP(stoic_model, 1);
fprintf('Conventional GEM growth (glc unlimited): %.4f 1/h\n', soln.f);

%% Stage 3 - Model Tuning: Sensitivity-based kcat tuning
% ecModel = loadEcModel('ecMMU1868a-pap_stage2.yml');
ecModel = setParam(ecModel, 'lb',  params.bioRxn, 0);
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);

ecModel = setProtPoolSize(ecModel, params.Ptot, params.f, 0.3);
[ecModel, tunedKcats] = sensitivityTuning(ecModel);

ecModel = setProtPoolSize(ecModel, params.Ptot, params.f, params.sigma);
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);
sol = solveLP(ecModel, 1);
fprintf('Stage 3 (post tuning) growth: %.4f 1/h\n', sol.f);

% Save the tuned model and the tuning report
saveEcModel(ecModel, 'ecMMU1868a-pap_stage3_tuned.yml');
writetable(struct2table(tunedKcats), ...
    fullfile(params.path, 'output', 'tunedKcats.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');

%% Stage 4 - Proteomics integration
protData = loadProtData(3);
ecModel  = fillEnzConcs(ecModel, protData);
ecModel  = constrainEnzConcs(ecModel);

f = calculateFfactor(ecModel, protData);
fprintf('Recalculated f factor (Sharma proteomics): %.4f\n', f);
%Calculated f factor: 0.0278

ecModel = setProtPoolSize(ecModel, params.Ptot, f, params.sigma);
ecModel = setParam(ecModel, 'lb',  params.bioRxn, 0);
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);
sol = solveLP(ecModel, 1);
fprintf('Stage 4 growth after proteomics: %.4f 1/h\n', sol.f);
saveEcModel(ecModel, 'ecMMU1868a-pap_stage4_prot.yml');

% Enzyme concentrations are flexibilized (increased) to recover max growth rate (gR_exp).
% Reset prot pool exchange to a safe value before flexibilization.
ecModel = setParam(ecModel, 'lb', 'prot_pool_exchange', -1000);
ecModel = setParam(ecModel, 'ub', 'prot_pool_exchange', 0);

[ecModel, flexEnz] = flexibilizeEnzConcs(ecModel, params.gR_exp, 10);

% Check growth
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);
solFlex = solveLP(ecModel, 1);
fprintf('Stage 4 growth after flexibilization: %.4f 1/h (target %.4f)\n', ...
    solFlex.f, params.gR_exp);

saveEcModel(ecModel, 'ecMMU1868a-pap_final.yml');

%% Stage 5 - Validation, enzyme usage, and ecModel vs conventional growth check
% Load the final ecModel and check maximum growth rate
ecModel = loadEcModel('ecMMU1868a-pap_final.yml');
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);
sol = solveLP(ecModel, 1);
ppIdx = getIndexes(ecModel, 'prot_pool_exchange', 'rxns');
fprintf('Final growth: %.4f 1/h\n', sol.f);
fprintf('prot_pool_exchange flux: %.6f  (lb=%.4g, ub=%.4g)\n', ...
    sol.x(ppIdx), ecModel.lb(ppIdx), ecModel.ub(ppIdx));

% Minimum protein at 99% of max growth
ecModel = setParam(ecModel, 'lb',  params.bioRxn, 0.99 * sol.f);
ecModel = setParam(ecModel, 'obj', 'prot_pool_exchange', 1);
solMinProt = solveLP(ecModel, 1);
fprintf('Min prot pool flux at 99%% growth: %.4f\n', solMinProt.x(ppIdx));

% Reset and report top enzyme usage
ecModel = setParam(ecModel, 'lb',  params.bioRxn, 0);
ecModel = setParam(ecModel, 'obj', params.bioRxn, 1);
sol = solveLP(ecModel, 1);
usageData   = enzymeUsage(ecModel, sol.x);
usageReport = reportEnzymeUsage(ecModel, usageData);
disp(usageReport.topAbsUsage);

% Compare ecModel vs conventional GEM
stoic_model = loadConventionalGEM();
stoic_model = setParam(stoic_model, 'obj', params.bioRxn, 1);
soln = solveLP(stoic_model, 1);
fprintf('Conventional growth: %.4f\n', soln.f);
fprintf('ecModel growth : %.4f\n', sol.f);

%% Save final ecModel as MAT for downstream SteadyCom (step 07_steadyCom_community)
ecModel = loadEcModel('ecMMU1868a-pap_final.yml');
save(fullfile(params.path, 'output', 'ecMMU1868a-pap_final.mat'), ...
    'ecModel', '-v7.3');
