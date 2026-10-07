% ============================================================================
% SteadyCom script for combining the three expanded enzyme-constrained cell
% models (ecMMU1868a-pap, ecMMU1868n-pre, and ecMMU1868n-post) into the
% enzyme-constrained tripartite synapse model (ecMMU1868-TPS).
%
% The three cell models are merged into a single SteadyCom community model.
% The extracellular compartment [e] is merged into the shared community pool
% [u] using the createMultipleSpeciesModel function. The synaptic cleft [syn]
% is merged into a shared synaptic pool [q] through reversible transfer
% reactions, enabling metabolite exchange across the synaptic cleft. Each
% species retains its own protein pool constraint (AST_, PRE_, and POST_)
% throughout community assembly.
%
% Equal-growth coupling, the defined medium, the protein-pool bound and the
% FBA objective are applied in this same script, so no intermediate
% ecTPS_build.mat is written.
%
% Note: createMultipleSpeciesModel drops the model.ec metadata field. The
% enzyme constraints survive in the S-matrix as 1/kcat coefficients, which is
% what the downstream analyses use.
%
% Output: ecMMU1868-TPS.mat
% ============================================================================

clc; clear;

%% Toolbox setup
COBRA_PATH  = 'C:/cobra';
RAVEN_PATH  = 'C:/RAVEN';
GECKO_PATH  = 'C:/GECKO';
SOLVER_PATH = 'C:/gurobi1202/win64/matlab';

restoredefaultpath; rehash toolboxcache;
set(0, 'DefaultFigureVisible', 'off');
addpath(genpath(COBRA_PATH));
addpath(genpath(RAVEN_PATH));
addpath(genpath(GECKO_PATH));
addpath(SOLVER_PATH);

initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');

%% 1) Paths
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

ecRoot   = fullfile(scriptDir, '..', '..', '06_enzyme_constraints_GECKO');
repoRoot = fullfile(scriptDir, '..', '..', '..');
outDir   = fullfile(repoRoot, 'models', '06_community_TPS', 'enzyme_constrained_tps');
if ~exist(outDir, 'dir'), mkdir(outDir); end

astFile  = fullfile(ecRoot, '1_ecMMU1868a-pap',  'output', 'ecMMU1868a-pap_final.mat');
preFile  = fullfile(ecRoot, '2_ecMMU1868n-pre',  'output', 'ecMMU1868n-pre_final.mat');
postFile = fullfile(ecRoot, '3_ecMMU1868n-post', 'output', 'ecMMU1868n-post_final.mat');
dstFile  = fullfile(outDir, 'ecMMU1868-TPS.mat');
for f = {astFile, preFile, postFile}
    assert(exist(f{1}, 'file') == 2, 'Missing %s', f{1});
end

%% 2) Load
fprintf('Loading ec models...\n');
matStruct = load(astFile);  AST  = pickEcModel(matStruct);
matStruct = load(preFile);  PRE  = pickEcModel(matStruct);
matStruct = load(postFile); POST = pickEcModel(matStruct);

AST = ensureBracketComp(AST); PRE = ensureBracketComp(PRE); POST = ensureBracketComp(POST);

fprintf('  Astrocyte : rxns=%d mets=%d genes=%d\n', numel(AST.rxns),  numel(AST.mets),  numel(AST.genes));
fprintf('  Pre-neuron: rxns=%d mets=%d genes=%d\n', numel(PRE.rxns),  numel(PRE.mets),  numel(PRE.genes));
fprintf('  Post-neuron:rxns=%d mets=%d genes=%d\n', numel(POST.rxns), numel(POST.mets), numel(POST.genes));

% confirm the new reactions are present BEFORE merging
newR = {'C04717tbulk_EXP_1','C04717tbulk_EXP_2','C04717tcFABP5_EXP_1','C04717tcFABP5_EXP_2'};
fprintf('\nNew reactions in the astrocyte ec model:\n');
for k = 1:numel(newR)
    fprintf('  %-24s %d\n', newR{k}, any(strcmp(AST.rxns, newR{k})));
end

%% 3) Clean empty cell fields
fieldsToClean = {'metFormulas','genes','grRules','metNames','rxnNames','subSystems','rules','geneNames'};
for f = fieldsToClean
    ff = f{1};
    if isfield(AST,ff),  AST.(ff)(cellfun(@isempty,AST.(ff)))   = {''}; end
    if isfield(PRE,ff),  PRE.(ff)(cellfun(@isempty,PRE.(ff)))   = {''}; end
    if isfield(POST,ff), POST.(ff)(cellfun(@isempty,POST.(ff))) = {''}; end
end

%% 4) Build community
desiredTags = {'AST_'; 'PRE_'; 'POST_'};
modArr = {AST, PRE, POST};
biomassNames = cell(3,1);
for k = 1:3
    bIdx = find(contains(modArr{k}.rxns, 'BIOMASS', 'IgnoreCase', true), 1);
    if ~isempty(bIdx), biomassNames{k} = modArr{k}.rxns{bIdx}; else, biomassNames{k} = ''; end
end
fprintf('\nBiomass reactions: %s | %s | %s\n', biomassNames{:});

modelCom = createMultipleSpeciesModel({AST; PRE; POST}, biomassNames, 'nameTagsModels', desiredTags);
modelCom.csense = repmat('E', 1, numel(modelCom.mets));
detectedTags = detectSpeciesPrefixes(modelCom.rxns, 3);
[modelCom.infoCom, modelCom.indCom] = getMultiSpeciesModelId(modelCom, detectedTags);
fprintf('Community: rxns=%d mets=%d\n', numel(modelCom.rxns), numel(modelCom.mets));

modelCom = trimSpuriousEXspecies(modelCom);
if ~isfield(modelCom,'genes') || isempty(modelCom.genes) || ...
   ~isfield(modelCom,'rxnGeneMat') || size(modelCom.rxnGeneMat,2)==0
    modelCom = restoreGeneFieldsFromSpecies(modelCom, {AST, PRE, POST}, detectedTags);
    fprintf('Rebuilt gene fields: total genes = %d\n', numel(modelCom.genes));
end

%% 5) Map species exchange bounds into [u]
mapToU = @(x) strrep(x, '[e]', '[u]');
getExchMet = @(mm) find(~cellfun('isempty', regexp(mm.mets, '\[e\]')));
metEx_AST = getExchMet(AST);  rxn_AST  = mapE_MetToExRxn(AST,  metEx_AST);
metEx_PRE = getExchMet(PRE);  rxn_PRE  = mapE_MetToExRxn(PRE,  metEx_PRE);
metEx_POST= getExchMet(POST); rxn_POST = mapE_MetToExRxn(POST, metEx_POST);

lb_AST = nan(numel(metEx_AST),1);  lb_PRE = nan(numel(metEx_PRE),1);  lb_POST = nan(numel(metEx_POST),1);
lb_AST(~isnan(rxn_AST))  = AST.lb(rxn_AST(~isnan(rxn_AST)));
lb_PRE(~isnan(rxn_PRE))  = PRE.lb(rxn_PRE(~isnan(rxn_PRE)));
lb_POST(~isnan(rxn_POST))= POST.lb(rxn_POST(~isnan(rxn_POST)));

[fA, iA_] = ismember(mapToU(AST.mets(metEx_AST)),  modelCom.infoCom.Mcom);
[fP, iP_] = ismember(mapToU(PRE.mets(metEx_PRE)),  modelCom.infoCom.Mcom);
[fQ, iQ_] = ismember(mapToU(POST.mets(metEx_POST)),modelCom.infoCom.Mcom);

modelCom.lb(modelCom.indCom.EXcom(:,1)) = 0;
modelCom.ub(modelCom.indCom.EXcom(:,1)) = 1e5;
for t = find(fA).', j = modelCom.indCom.EXcom(iA_(t),1); modelCom.lb(j) = min(modelCom.lb(j), lb_AST(t)); end
for t = find(fP).', j = modelCom.indCom.EXcom(iP_(t),1); modelCom.lb(j) = min(modelCom.lb(j), lb_PRE(t)); end
for t = find(fQ).', j = modelCom.indCom.EXcom(iQ_(t),1); modelCom.lb(j) = min(modelCom.lb(j), lb_POST(t)); end

[modelCom.infoCom, modelCom.indCom] = getMultiSpeciesModelId(modelCom, detectedTags);
modelCom = trimSpuriousEXspecies(modelCom);

%% 6) Shared synaptic cleft pool [q]
modelCom = addSharedPoolForCompartment(modelCom, detectedTags, 'syn', 'q');

%% 7) Equal-flux biomass coupling
iA = find(strcmp(modelCom.rxns, 'AST_BIOMASS_reaction'),  1);
iP = find(strcmp(modelCom.rxns, 'PRE_BIOMASS_reaction'),  1);
iQ = find(strcmp(modelCom.rxns, 'POST_BIOMASS_reaction'), 1);
assert(~isempty(iA) && ~isempty(iP) && ~isempty(iQ), 'biomass reactions not found');

nRxn = numel(modelCom.rxns);
modelCom.S = [modelCom.S; sparse(1,[iA iP],[1 -1],1,nRxn); sparse(1,[iA iQ],[1 -1],1,nRxn)];
modelCom.b = [modelCom.b; 0; 0];
if isfield(modelCom,'csense'), modelCom.csense = [modelCom.csense(:); 'E'; 'E']; end
modelCom.mets = [modelCom.mets; {'couple_AST_PRE'; 'couple_AST_POST'}];
modelCom = padMetFields(modelCom, 2);
fprintf('Added equal-flux coupling.\n');

%% 8) Defined minimal medium
medium = { ...
  'EX_h2o[u]' -1000 1000; 'EX_h[u]' -1000 1000; 'EX_na1[u]' -1000 1000;
  'EX_k[u]' -1000 1000;   'EX_cl[u]' -1000 1000; 'EX_ca2[u]' -1000 1000;
  'EX_nh4[u]' -10 1000;   'EX_pi[u]' -10 1000;   'EX_so4[u]' -10 1000;
  'EX_o2[u]' -10 1000;    'EX_co2[u]' 0 1000;    'EX_glc__D[u]' -10 1000;
  'EX_gln__L[u]' 0 1000;  'EX_lac__L[u]' 0 1000; 'EX_hdca[u]' 0 1000;
  'EX_his__L[u]' -0.186207 1000; 'EX_ile__L[u]' -0.186207 1000;
  'EX_leu__L[u]' -0.186207 1000; 'EX_lys__L[u]' -0.186207 1000;
  'EX_met__L[u]' -0.155172 1000; 'EX_phe__L[u]' -0.279310 1000;
  'EX_thr__L[u]' -0.124138 1000; 'EX_trp__L[u]' -0.034483 1000;
  'EX_val__L[u]' -0.155172 1000 };

exIdx = find(~cellfun('isempty', regexp(modelCom.rxns, '^EX_.*\[u\]$')));
modelCom.lb(exIdx) = 0; modelCom.ub(exIdx) = 1000;
nOpened = 0;
for i = 1:size(medium,1)
    j = find(strcmp(modelCom.rxns, medium{i,1}), 1);
    if ~isempty(j), modelCom.lb(j) = medium{i,2}; modelCom.ub(j) = medium{i,3}; nOpened = nOpened + 1; end
end
fprintf('Medium: %d/%d exchanges opened.\n', nOpened, size(medium,1));

%% 9) prot_pool bound
for tag = {'AST','PRE','POST'}
    j = find(strcmp(modelCom.rxns, [tag{1} '_prot_pool_exchange']), 1);
    if ~isempty(j), modelCom.lb(j) = -10; end
end
fprintf('prot_pool_exchange lb = -10 for all three species.\n');

%% 10) Objective
modelCom.c = zeros(numel(modelCom.rxns),1);
modelCom.c(iA) = 1;
modelCom.osenseStr = 'max';

%% 11) VERIFY the new FABP reactions survived the merge
fprintf('\n=========================================================\n');
fprintf('VERIFICATION -- new FABP reactions in the community model\n');
fprintf('=========================================================\n');
checkR = {'AST_C04717tbulk_EXP_1','AST_C04717tbulk_EXP_2', ...
          'AST_C04717tbulk_REV_EXP_1','AST_C04717tbulk_REV_EXP_2', ...
          'AST_C04717tcFABP5_EXP_1','AST_C04717tcFABP5_EXP_2', ...
          'AST_C04717tc','AST_C04717td','AST_C04717td_REV'};
for k = 1:numel(checkR)
    j = find(strcmp(modelCom.rxns, checkR{k}), 1);
    if isempty(j)
        fprintf('  %-30s MISSING\n', checkR{k});
    else
        protRows = find(modelCom.S(:,j) ~= 0 & contains(modelCom.mets, 'prot_'));
        if isempty(protRows)
            fprintf('  %-30s present  lb=%g ub=%g  (no prot cost)\n', checkR{k}, modelCom.lb(j), modelCom.ub(j));
        else
            fprintf('  %-30s present  lb=%g ub=%g  cost=%.6f on %s\n', ...
                checkR{k}, modelCom.lb(j), modelCom.ub(j), ...
                full(-modelCom.S(protRows(1), j)), modelCom.mets{protRows(1)});
        end
    end
end
jb = find(strcmp(modelCom.mets, 'AST_C04717bulk[c]'), 1);
if isempty(jb)
    fprintf('  AST_C04717bulk[c] MISSING\n');
else
    fprintf('  AST_C04717bulk[c] present, in %d reactions\n', nnz(modelCom.S(jb,:)));
end

%% 12) Save, then verify growth
if exist(dstFile,'file'), fileattrib(dstFile,'+w','a'); delete(dstFile); end
save(dstFile, 'modelCom', '-v7.3');
fprintf('\nSaved %s\n', dstFile);

fprintf('\nOptimising ...\n');
sol = optimizeCbModel(modelCom, 'max');
fprintf('  status = %d\n', sol.stat);
fprintf('  mu_max = %.6f 1/h\n', sol.f);
fprintf('  AST_BM = %.6f\n  PRE_BM = %.6f\n  POST_BM = %.6f\n', ...
    sol.x(iA), sol.x(iP), sol.x(iQ));
fprintf('\nCommunity: rxns=%d mets=%d genes=%d\n', ...
    numel(modelCom.rxns), numel(modelCom.mets), numel(modelCom.genes));


%% ===================== local functions (from Tripartite_ecCom.m) =========
function ecModel = pickEcModel(matStruct)
    if isfield(matStruct,'ecModel'), ecModel = matStruct.ecModel;
    elseif isfield(matStruct,'model'), ecModel = matStruct.model;
    elseif isfield(matStruct,'modelEC'), ecModel = matStruct.modelEC;
    else
        fields = fieldnames(matStruct);
        for fi = 1:numel(fields)
            val = matStruct.(fields{fi});
            if isstruct(val) && isfield(val,'S') && isfield(val,'rxns'), ecModel = val; return; end
        end
        error('No COBRA-like struct found in MAT file');
    end
end

function model = ensureBracketComp(model)
    if any(contains(model.mets,'[')), return; end
    if ~isfield(model,'comps') || isempty(model.comps), return; end
    newMets = model.mets;
    for mi = 1:numel(newMets)
        metID = newMets{mi};
        for ci = 1:numel(model.comps)
            suffix = ['_', model.comps{ci}];
            if endsWith(metID, suffix)
                newMets{mi} = [metID(1:end-numel(suffix)), '[', model.comps{ci}, ']']; break;
            end
        end
    end
    model.mets = newMets;
end

function modelCom = trimSpuriousEXspecies(modelCom)
    if ~isfield(modelCom,'infoCom') || ~isfield(modelCom.infoCom,'spAbbr'), return; end
    keep = ~strcmp(modelCom.infoCom.spAbbr,'EX_');
    if any(~keep)
        modelCom.infoCom.spAbbr = modelCom.infoCom.spAbbr(keep);
        if isfield(modelCom.infoCom,'spName'), modelCom.infoCom.spName = modelCom.infoCom.spName(keep); end
        nSp = numel(modelCom.infoCom.spAbbr);
        if size(modelCom.indCom.EXsp,2) > nSp, modelCom.indCom.EXsp = modelCom.indCom.EXsp(:,1:nSp); end
    end
end

function rxnIdx = mapE_MetToExRxn(model, extracellularMets)
    singleMetRxns = find(sum(model.S ~= 0, 1) == 1);
    [metRows, ~] = find(model.S(:, singleMetRxns));
    isExtracellular = false(size(singleMetRxns));
    for col = 1:numel(singleMetRxns)
        compTok = regexp(model.mets{metRows(col)}, '\[(.)\]', 'tokens', 'once');
        isExtracellular(col) = ~isempty(compTok) && strcmp(compTok{1},'e');
    end
    extracellularExRxns = singleMetRxns(isExtracellular);
    rxnIdx = nan(numel(extracellularMets),1);
    for mi = 1:numel(extracellularMets)
        matchCols = find(model.S(extracellularMets(mi), extracellularExRxns) ~= 0);
        if numel(matchCols) == 1, rxnIdx(mi) = extracellularExRxns(matchCols); end
    end
end

function tags = detectSpeciesPrefixes(rxnIDs, nSpecies)
    prefixTokens = regexp(rxnIDs, '^([A-Za-z0-9]+_)', 'tokens', 'once');
    prefixes = prefixTokens(~cellfun(@isempty, prefixTokens));
    prefixes = cellfun(@(tk) tk{1}, prefixes, 'UniformOutput', false);
    [~, firstIdx] = unique(prefixes, 'stable');
    prefixes = prefixes(sort(firstIdx));
    if isempty(prefixes) || numel(prefixes) < nSpecies
        prefixes = arrayfun(@(i) sprintf('model%d_',i), 1:nSpecies, 'UniformOutput', false);
    end
    tags = prefixes(:);
end

function modelCom = addSharedPoolForCompartment(modelCom, tags, compFrom, compPool)
    fromPattern = ['\[', compFrom, '\]'];
    addedPoolMets = containers.Map;
    for si = 1:numel(tags)
        tag = tags{si};
        speciesMetIdx = find(startsWith(modelCom.mets, tag) & ...
                        ~cellfun('isempty', regexp(modelCom.mets, fromPattern)));
        for mi = 1:numel(speciesMetIdx)
            speciesMet = modelCom.mets{speciesMetIdx(mi)};
            poolMet = regexprep(erase(speciesMet, tag), fromPattern, ['[', compPool, ']']);
            if ~isKey(addedPoolMets, poolMet)
                addedPoolMets(poolMet) = 1;
                if findMetIDs(modelCom, poolMet) == 0, modelCom = addMetabolite(modelCom, poolMet); end
                exRxn = ['EXC_', regexprep(poolMet, ['\[',compPool,'\]'], ''), '[', compPool, ']'];
                if findRxnIDs(modelCom, exRxn) == 0
                    modelCom = addReaction(modelCom, exRxn, 'metaboliteList', {poolMet}, ...
                        'stoichCoeffList', -1, 'lowerBound', 0, 'upperBound', 0);
                end
            end
            transferRxn = ['ICLEFT_', tag, regexprep(poolMet, ['\[',compPool,'\]'], ''), ...
                           '_', compFrom, '_to_', compPool];
            if findRxnIDs(modelCom, transferRxn) == 0
                modelCom = addReaction(modelCom, transferRxn, ...
                    'metaboliteList', {speciesMet, poolMet}, 'stoichCoeffList', [1,-1], ...
                    'lowerBound', -1e5, 'upperBound', 1e5);
            end
        end
    end
    fprintf('Added shared cleft pool [%s] from [%s]\n', compPool, compFrom);
end

function modelCom = restoreGeneFieldsFromSpecies(modelCom, speciesModels, tags)
    nRxnsCom = numel(modelCom.rxns); nSpecies = numel(speciesModels);
    nGenesPerSp = zeros(nSpecies,1);
    for si = 1:nSpecies
        spModel = speciesModels{si};
        if isfield(spModel,'genes') && ~isempty(spModel.genes), nGenesPerSp(si) = numel(spModel.genes); end
    end
    nGenesTotal = sum(nGenesPerSp);
    allGenes = cell(nGenesTotal,1); allGeneNames = cell(nGenesTotal,1);
    geneColStart = [0; cumsum(nGenesPerSp(1:end-1))] + 1;
    rxnGeneMatCom = sparse(nRxnsCom, nGenesTotal);
    rules = repmat({''}, nRxnsCom, 1); grRules = repmat({''}, nRxnsCom, 1);
    comRxnMap = containers.Map(modelCom.rxns, num2cell(1:nRxnsCom));
    for si = 1:nSpecies
        spModel = speciesModels{si}; tag = tags{si};
        gStart = geneColStart(si); nGenes = nGenesPerSp(si);
        for gi = 1:nGenes
            allGenes{gStart+gi-1} = [tag, spModel.genes{gi}];
            if isfield(spModel,'geneNames') && numel(spModel.geneNames) >= gi && ~isempty(spModel.geneNames{gi})
                allGeneNames{gStart+gi-1} = spModel.geneNames{gi};
            else
                allGeneNames{gStart+gi-1} = spModel.genes{gi};
            end
        end
        if nGenes > 0
            escapedGenes = cellfun(@regexpescape, spModel.genes, 'UniformOutput', false);
            [~, lenOrder] = sort(cellfun(@numel, spModel.genes), 'descend');
            genePattern = ['(?<![A-Za-z0-9_])(?:', strjoin(escapedGenes(lenOrder),'|'), ')(?![A-Za-z0-9_])'];
        else
            genePattern = '';
        end
        geneIdxOffset = gStart - 1;
        for ri = 1:numel(spModel.rxns)
            comRxnID = [tag, spModel.rxns{ri}];
            if ~isKey(comRxnMap, comRxnID), continue; end
            comRxnIdx = comRxnMap(comRxnID);
            if nGenes > 0 && isfield(spModel,'rxnGeneMat') && ~isempty(spModel.rxnGeneMat)
                geneCols = find(spModel.rxnGeneMat(ri,:));
                if ~isempty(geneCols), rxnGeneMatCom(comRxnIdx, gStart-1+geneCols) = spModel.rxnGeneMat(ri, geneCols); end
            end
            if isfield(spModel,'rules') && numel(spModel.rules) >= ri && ~isempty(spModel.rules{ri})
                rules{comRxnIdx} = remapRuleIndices(spModel.rules{ri}, geneIdxOffset);
            end
            if isfield(spModel,'grRules') && numel(spModel.grRules) >= ri && ~isempty(spModel.grRules{ri})
                if ~isempty(genePattern)
                    grRules{comRxnIdx} = regexprep(spModel.grRules{ri}, genePattern, [tag, '$0']);
                else
                    grRules{comRxnIdx} = spModel.grRules{ri};
                end
            end
        end
        fprintf('  Gene mapping for %s complete (%d rxns, %d genes)\n', tag, numel(spModel.rxns), nGenes);
    end
    modelCom.genes = allGenes; modelCom.geneNames = allGeneNames;
    modelCom.rxnGeneMat = rxnGeneMatCom;
    if any(~cellfun(@isempty, rules)),   modelCom.rules   = rules;   end
    if any(~cellfun(@isempty, grRules)), modelCom.grRules = grRules; end
end

function out = remapRuleIndices(ruleIn, geneIdxOffset)
    if geneIdxOffset == 0, out = ruleIn; return; end
    idxTokens = regexp(ruleIn, 'x\((\d+)\)', 'tokens');
    tokStart  = regexp(ruleIn, 'x\((\d+)\)', 'start');
    tokEnd    = regexp(ruleIn, 'x\((\d+)\)', 'end');
    if isempty(idxTokens), out = ruleIn; return; end
    out = ''; prevEnd = 0;
    for k = 1:numel(idxTokens)
        out = [out, ruleIn(prevEnd+1:tokStart(k)-1), sprintf('x(%d)', str2double(idxTokens{k}{1}) + geneIdxOffset)];
        prevEnd = tokEnd(k);
    end
    out = [out, ruleIn(prevEnd+1:end)];
end

function str = regexpescape(str)
    metaChars = '.[{}()\^$|*+?-\]';
    for ch = metaChars, str = strrep(str, ch, ['\', ch]); end
end

function model = padMetFields(model, k)
    if isfield(model,'metNames'),   model.metNames   = [model.metNames;   repmat({''},k,1)]; end
    if isfield(model,'metFormulas'),model.metFormulas= [model.metFormulas;repmat({''},k,1)]; end
    if isfield(model,'metCharges'), model.metCharges = [model.metCharges; zeros(k,1)]; end
    if isfield(model,'metComps'),   model.metComps   = [model.metComps;   zeros(k,1)]; end
end
