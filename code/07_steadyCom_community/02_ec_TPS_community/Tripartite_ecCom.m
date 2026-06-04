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
% Output: ecTPS_build.mat
% ============================================================================

clc; clear;

cd(fileparts(which(mfilename)));

%% Toolbox setup
% Edit toolbox paths
% COBRA Toolbox : https://opencobra.github.io/cobratoolbox/
% RAVEN         : https://github.com/SysBioChalmers/RAVEN
% GECKO         : https://github.com/SysBioChalmers/GECKO
% Gurobi        : https://www.gurobi.com
COBRA_PATH = 'C:/cobra';
RAVEN_PATH = 'C:/RAVEN';
GECKO_PATH = 'C:/GECKO';
SOLVER_PATH = 'C:/gurobi1202/win64/matlab';

restoredefaultpath; rehash toolboxcache;
addpath(genpath(COBRA_PATH));
addpath(genpath(RAVEN_PATH));
addpath(genpath(GECKO_PATH));
addpath(SOLVER_PATH);

initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');

%% 1) Locate the three ecModel .mat files
% Find the .mat files archive in zenodo or convert them final .yml ecModels
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

ecRoot = fullfile(scriptDir, '..', '..', '06_enzyme_constraints_GECKO');
astFile = fullfile(ecRoot, '1_ecMMU1868a-pap', 'output', 'ecMMU1868a-pap_final.mat');
preFile = fullfile(ecRoot, '2_ecMMU1868n-pre', 'output', 'ecMMU1868n-pre_final.mat');
postFile = fullfile(ecRoot, '3_ecMMU1868n-post', 'output', 'ecMMU1868n-post_final.mat');

assert(exist(astFile, 'file') == 2, 'Missing %s', astFile);
assert(exist(preFile, 'file') == 2, 'Missing %s', preFile);
assert(exist(postFile, 'file') == 2, 'Missing %s', postFile);

repoRoot = fullfile(scriptDir, '..', '..', '..');
outDir   = fullfile(repoRoot, 'models', '06_community_TPS', 'enzyme_constrained_tps');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% 2) Load ecModel structs
fprintf('Loading ec models...\n');
matStruct = load(astFile); AST  = pickEcModel(matStruct);
matStruct = load(preFile); PRE  = pickEcModel(matStruct);
matStruct = load(postFile); POST = pickEcModel(matStruct);

% Ensure metabolite IDs use the [comp] suffix expected by SteadyCom helpers
AST = ensureBracketComp(AST);
PRE = ensureBracketComp(PRE);
POST = ensureBracketComp(POST);

fprintf('  Astrocyte (PAP):   rxns=%d, mets=%d, genes=%d\n', ...
    numel(AST.rxns), numel(AST.mets), numel(AST.genes));
fprintf('  Pre-neuron:        rxns=%d, mets=%d, genes=%d\n', ...
    numel(PRE.rxns), numel(PRE.mets), numel(PRE.genes));
fprintf('  Post-neuron:       rxns=%d, mets=%d, genes=%d\n', ...
    numel(POST.rxns), numel(POST.mets), numel(POST.genes));

%% 3) Clean empty cell fields (required by createMultipleSpeciesModel)
fieldsToClean = {'metFormulas','genes','grRules','metNames','rxnNames', ...
                 'subSystems','rules','geneNames'};
for f = fieldsToClean
    ff = f{1};
    if isfield(AST,ff),  AST.(ff)(cellfun(@isempty,AST.(ff)))   = {''}; end
    if isfield(PRE,ff),  PRE.(ff)(cellfun(@isempty,PRE.(ff)))   = {''}; end
    if isfield(POST,ff), POST.(ff)(cellfun(@isempty,POST.(ff))) = {''}; end
end

%% 4) Build multi-species community model
desiredTags = {'AST_'; 'PRE_'; 'POST_'};

biomassNames = cell(3,1);
modArr = {AST, PRE, POST};
for k = 1:3
    bIdx = find(contains(modArr{k}.rxns, 'BIOMASS', 'IgnoreCase', true), 1);
    if ~isempty(bIdx)
        biomassNames{k} = modArr{k}.rxns{bIdx};
    else
        biomassNames{k} = '';
    end
end
fprintf('\nBiomass reactions: %s | %s | %s\n', biomassNames{:});

modelCom = createMultipleSpeciesModel({AST; PRE; POST}, biomassNames, ...
    'nameTagsModels', desiredTags);
modelCom.csense = repmat('E', 1, numel(modelCom.mets));

detectedTags = detectSpeciesPrefixes(modelCom.rxns, 3);
[modelCom.infoCom, modelCom.indCom] = getMultiSpeciesModelId(modelCom, detectedTags);

fprintf('Detected rxn prefixes: %s, %s, %s\n', detectedTags{:});
fprintf('Community: rxns=%d, mets=%d\n', ...
    numel(modelCom.rxns), numel(modelCom.mets));

%% 5) Remove EX_ "species" if present
modelCom = trimSpuriousEXspecies(modelCom);
fprintf('Species count: %d\n', numel(modelCom.infoCom.spAbbr));

%% 6) Rebuild gene fields from individual species models
if ~isfield(modelCom,'genes') || isempty(modelCom.genes) || ...
   ~isfield(modelCom,'rxnGeneMat') || size(modelCom.rxnGeneMat,2)==0
    modelCom = restoreGeneFieldsFromSpecies(modelCom, {AST, PRE, POST}, detectedTags);
    fprintf('Rebuilt gene fields: total genes = %d\n', numel(modelCom.genes));
end

%% 7) Map species exchange bounds into community [u] exchanges
mapToU  = @(metIDs) strrep(metIDs, '[e]', '[u]');
getExchMet = @(model) find(~cellfun('isempty', regexp(model.mets, '\[e\]')));

metEx_AST = getExchMet(AST);   rxn_AST = mapE_MetToExRxn(AST,  metEx_AST);
metEx_PRE = getExchMet(PRE);   rxn_PRE = mapE_MetToExRxn(PRE,  metEx_PRE);
metEx_POST = getExchMet(POST);  rxn_POST = mapE_MetToExRxn(POST, metEx_POST);

lb_AST = nan(numel(metEx_AST),1);
lb_PRE = nan(numel(metEx_PRE),1);
lb_POST = nan(numel(metEx_POST),1);
lb_AST(~isnan(rxn_AST))  = AST.lb(rxn_AST(~isnan(rxn_AST)));
lb_PRE(~isnan(rxn_PRE)) = PRE.lb(rxn_PRE(~isnan(rxn_PRE)));
lb_POST(~isnan(rxn_POST)) = POST.lb(rxn_POST(~isnan(rxn_POST)));

[found_AST, idxComMet_AST] = ismember(mapToU(AST.mets(metEx_AST)), modelCom.infoCom.Mcom);
[found_PRE, idxComMet_PRE] = ismember(mapToU(PRE.mets(metEx_PRE)), modelCom.infoCom.Mcom);
[found_POST, idxComMet_POST]= ismember(mapToU(POST.mets(metEx_POST)), modelCom.infoCom.Mcom);

modelCom.lb(modelCom.indCom.EXcom(:,1)) = 0;
modelCom.ub(modelCom.indCom.EXcom(:,1)) = 1e5;
for t = find(found_AST).'
    j = modelCom.indCom.EXcom(idxComMet_AST(t),1);
    modelCom.lb(j) = min(modelCom.lb(j), lb_AST(t));
end
for t = find(found_PRE).'
    j = modelCom.indCom.EXcom(idxComMet_PRE(t),1);
    modelCom.lb(j) = min(modelCom.lb(j), lb_PRE(t));
end
for t = find(found_POST).'
    j = modelCom.indCom.EXcom(idxComMet_POST(t),1);
    modelCom.lb(j) = min(modelCom.lb(j), lb_POST(t));
end
fprintf('Community exchange bounds mapped from species models.\n');

%% 7b) Refresh infoCom/indCom and re-trim spurious EX_ species
[modelCom.infoCom, modelCom.indCom] = getMultiSpeciesModelId(modelCom, detectedTags);
modelCom = trimSpuriousEXspecies(modelCom);
fprintf('infoCom/indCom refreshed (EXcom: %d, EXsp: %d x %d)\n', ...
    size(modelCom.indCom.EXcom,1), size(modelCom.indCom.EXsp,1), size(modelCom.indCom.EXsp,2));

%% 8) Add the shared synaptic cleft pool [q]
modelCom = addSharedPoolForCompartment(modelCom, detectedTags, 'syn', 'q');

synPoolMets = modelCom.mets(contains(modelCom.mets, '[q]'));
fprintf('\nShared pool [q] metabolites: %d\n', numel(synPoolMets));

%% 9) Community model summary
fprintf('\n EC COMMUNITY MODEL SUMMARY: \n');
fprintf('Total reactions: %d\n', numel(modelCom.rxns));
fprintf('Total metabolites: %d\n', numel(modelCom.mets));
fprintf('Total genes: %d\n', numel(modelCom.genes));

fprintf('Species: %s\n', strjoin(modelCom.infoCom.spAbbr, ', '));
fprintf('EXcom: %d, EXsp: %d x %d\n', ...
    size(modelCom.indCom.EXcom,1), size(modelCom.indCom.EXsp,1), size(modelCom.indCom.EXsp,2));
icleft_rxns = modelCom.rxns(startsWith(modelCom.rxns, 'ICLEFT_'));
fprintf('Synaptic coupling transfers:   %d\n', numel(icleft_rxns));

%% 10) Save community model
matFile = fullfile(outDir, 'ecTPS_build.mat');
if exist(matFile, 'file'); delete(matFile); end
save(matFile, 'modelCom', '-v7.3');
fprintf('\nSaved %s\n', matFile);
 
%% Local Functions
% Required by the script to run properly. 

function ecModel = pickEcModel(matStruct)
% Return the COBRA-like struct stored in a loaded ecModel .mat file.
    if isfield(matStruct, 'ecModel'), ecModel = matStruct.ecModel;
    elseif isfield(matStruct, 'model'), ecModel = matStruct.model;
    elseif isfield(matStruct, 'modelEC'), ecModel = matStruct.modelEC;
    else
        fields = fieldnames(matStruct);
        for fi = 1:numel(fields)
            val = matStruct.(fields{fi});
            if isstruct(val) && isfield(val, 'S') && isfield(val, 'rxns')
                ecModel = val; return;
            end
        end
        error('No COBRA-like struct found in MAT file');
    end
end


function model = ensureBracketComp(model)
% Convert metabolite IDs metName_comp -> metName[comp] when needed.
% GECKO ecModels keep underscore-compartment naming; SteadyCom expects [comp].
    if any(contains(model.mets, '['))
        return;   % already bracketed
    end
    if ~isfield(model, 'comps') || isempty(model.comps)
        return;
    end
    newMets = model.mets;
    for mi = 1:numel(newMets)
        metID = newMets{mi};
        for ci = 1:numel(model.comps)
            suffix = ['_', model.comps{ci}];
            if endsWith(metID, suffix)
                newMets{mi} = [metID(1:end-numel(suffix)), '[', model.comps{ci}, ']'];
                break;
            end
        end
    end
    model.mets = newMets;
end


function modelCom = trimSpuriousEXspecies(modelCom)
% Remove a unintended 'EX_' entry that COBRA sometimes lists as a species.
    if ~isfield(modelCom, 'infoCom') || ~isfield(modelCom.infoCom, 'spAbbr')
        return;
    end
    keep = ~strcmp(modelCom.infoCom.spAbbr, 'EX_');
    if any(~keep)
        modelCom.infoCom.spAbbr = modelCom.infoCom.spAbbr(keep);
        if isfield(modelCom.infoCom, 'spName')
            modelCom.infoCom.spName = modelCom.infoCom.spName(keep);
        end
        nSp = numel(modelCom.infoCom.spAbbr);
        if size(modelCom.indCom.EXsp, 2) > nSp
            modelCom.indCom.EXsp = modelCom.indCom.EXsp(:, 1:nSp);
        end
    end
end


function rxnIdx = mapE_MetToExRxn(model, extracellularMets)
% For each extracellular metabolite, find its exchange reaction (NaN if none).
    singleMetRxns = find(sum(model.S ~= 0, 1) == 1);

    [metRows, ~] = find(model.S(:, singleMetRxns));
    isExtracellular = false(size(singleMetRxns));
    for col = 1:numel(singleMetRxns)
        metRow  = metRows(col);
        compTok = regexp(model.mets{metRow}, '\[(.)\]', 'tokens', 'once');
        isExtracellular(col) = ~isempty(compTok) && strcmp(compTok{1}, 'e');
    end
    extracellularExRxns = singleMetRxns(isExtracellular);

    rxnIdx = nan(numel(extracellularMets), 1);
    for mi = 1:numel(extracellularMets)
        metRow    = extracellularMets(mi);
        matchCols = find(model.S(metRow, extracellularExRxns) ~= 0);
        if numel(matchCols) == 1
            rxnIdx(mi) = extracellularExRxns(matchCols);
        end
    end
end


function tags = detectSpeciesPrefixes(rxnIDs, nSpecies)
% Extract the species prefixes COBRA assigned (e.g. AST_, PRE_, POST_).
    prefixTokens = regexp(rxnIDs, '^([A-Za-z0-9]+_)', 'tokens', 'once');
    prefixes = prefixTokens(~cellfun(@isempty, prefixTokens));
    prefixes = cellfun(@(tk) tk{1}, prefixes, 'UniformOutput', false);
    [~, firstIdx] = unique(prefixes, 'stable');
    prefixes = prefixes(sort(firstIdx));
    if isempty(prefixes) || numel(prefixes) < nSpecies
        prefixes = arrayfun(@(idx) sprintf('model%d_', idx), 1:nSpecies, ...
            'UniformOutput', false);
    end
    tags = prefixes(:);
end


function modelCom = addSharedPoolForCompartment(modelCom, tags, compFrom, compPool)
% Create a shared metabolite pool [compPool] for an inter-cellular
% compartment [compFrom], linking each species' metabolites to the pool
% via a reversible coupling transfer (e.g. AST_met[syn] <-> met[q]).
    fromPattern = ['\[', compFrom, '\]'];
    addedPoolMets = containers.Map;
    for si = 1:numel(tags)
        tag = tags{si};
        speciesMetIdx = find(startsWith(modelCom.mets, tag) & ...
                        ~cellfun('isempty', regexp(modelCom.mets, fromPattern)));
        for mi = 1:numel(speciesMetIdx)
            speciesMet = modelCom.mets{speciesMetIdx(mi)};
            metNoTag = erase(speciesMet, tag);
            poolMet  = regexprep(metNoTag, fromPattern, ['[', compPool, ']']);
            if ~isKey(addedPoolMets, poolMet)
                addedPoolMets(poolMet) = 1;
                if findMetIDs(modelCom, poolMet) == 0
                    modelCom = addMetabolite(modelCom, poolMet);
                end
                exRxn = ['EXC_', regexprep(poolMet, ['\[', compPool, '\]'], ''), ...
                         '[', compPool, ']'];
                if findRxnIDs(modelCom, exRxn) == 0
                    modelCom = addReaction(modelCom, exRxn, ...
                        'metaboliteList', {poolMet}, ...
                        'stoichCoeffList', -1, ...
                        'lowerBound', 0, 'upperBound', 0);
                end
            end
            transferRxn = ['ICLEFT_', tag, ...
                     regexprep(poolMet, ['\[', compPool, '\]'], ''), ...
                     '_', compFrom, '_to_', compPool];
            if findRxnIDs(modelCom, transferRxn) == 0
                modelCom = addReaction(modelCom, transferRxn, ...
                    'metaboliteList', {speciesMet, poolMet}, ...
                    'stoichCoeffList', [1, -1], ...
                    'lowerBound', -1e5, 'upperBound', 1e5);
            end
        end
    end
    fprintf('Added shared cleft pool [%s] from species compartment [%s]\n', ...
        compPool, compFrom);
end


function modelCom = restoreGeneFieldsFromSpecies(modelCom, speciesModels, tags)
% Rebuild genes, rxnGeneMat, rules, and grRules (with species prefixes)
% from the individual models, since createMultipleSpeciesModel can drop them.
    nRxnsCom = numel(modelCom.rxns);
    nSpecies = numel(speciesModels);
    nGenesPerSp = zeros(nSpecies, 1);
    for si = 1:nSpecies
        spModel = speciesModels{si};
        if isfield(spModel,'genes') && ~isempty(spModel.genes)
            nGenesPerSp(si) = numel(spModel.genes);
        else
            nGenesPerSp(si) = 0;
        end
    end
    nGenesTotal = sum(nGenesPerSp);
    allGenes  = cell(nGenesTotal, 1);
    allGeneNames = cell(nGenesTotal, 1);
    geneColStart = [0; cumsum(nGenesPerSp(1:end-1))] + 1;
    rxnGeneMatCom = sparse(nRxnsCom, nGenesTotal);
    rules = repmat({''}, nRxnsCom, 1);
    grRules = repmat({''}, nRxnsCom, 1);
    comRxnMap = containers.Map(modelCom.rxns, num2cell(1:nRxnsCom));
    for si = 1:nSpecies
        spModel = speciesModels{si};  tag = tags{si};
        gStart  = geneColStart(si);   nGenes = nGenesPerSp(si);
        for gi = 1:nGenes
            allGenes{gStart+gi-1} = [tag, spModel.genes{gi}];
            if isfield(spModel,'geneNames') && numel(spModel.geneNames) >= gi && ~isempty(spModel.geneNames{gi})
                allGeneNames{gStart+gi-1} = spModel.geneNames{gi};
            else
                allGeneNames{gStart+gi-1} = spModel.genes{gi};
            end
        end
        % One combined regex to prefix all genes at once.
        if nGenes > 0
            escapedGenes = cellfun(@regexpescape, spModel.genes, 'UniformOutput', false);
            [~, lenOrder] = sort(cellfun(@numel, spModel.genes), 'descend');
            genePattern = ['(?<![A-Za-z0-9_])(?:', strjoin(escapedGenes(lenOrder), '|'), ')(?![A-Za-z0-9_])'];
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
        fprintf('  Gene mapping for %s complete (%d rxns, %d genes)\n', ...
            tag, numel(spModel.rxns), nGenes);
    end
    modelCom.genes = allGenes;
    modelCom.geneNames = allGeneNames;
    modelCom.rxnGeneMat = rxnGeneMatCom;
    if any(~cellfun(@isempty, rules)),   modelCom.rules   = rules;   end
    if any(~cellfun(@isempty, grRules)), modelCom.grRules = grRules; end
end


function out = remapRuleIndices(ruleIn, offset)
% Shift each gene index x(N) to x(N+offset) for community indexing.
    if offset == 0, out = ruleIn; return; end
    idxTokens = regexp(ruleIn, 'x\((\d+)\)', 'tokens');
    tokStart = regexp(ruleIn, 'x\((\d+)\)', 'start');
    tokEnd = regexp(ruleIn, 'x\((\d+)\)', 'end');
    if isempty(idxTokens), out = ruleIn; return; end
    out = ''; prevEnd = 0;
    for k = 1:numel(idxTokens)
        out = [out, ruleIn(prevEnd+1:tokStart(k)-1), ...
               sprintf('x(%d)', str2double(idxTokens{k}{1}) + offset)];
        prevEnd = tokEnd(k);
    end
    out = [out, ruleIn(prevEnd+1:end)];
end


function str = regexpescape(str)
% Escape regex metacharacters in a string.
    metaChars = '.[{}()\^$|*+?-\]';
    for ch = metaChars, str = strrep(str, ch, ['\', ch]); end
end
