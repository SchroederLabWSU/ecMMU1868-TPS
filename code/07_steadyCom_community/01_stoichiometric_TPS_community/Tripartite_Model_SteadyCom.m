% ============================================================================
% SteadyCom script for combining the three expanded stoichiometric cell models
% (iMMU1868a-pap, iMMU1868n-pre and iMMU1868n-post) into stoichiometric 
% tripartite synapse model (iMMU1868-TPS)
%
% The three cell models are merged into a single SteadyCom community model. The
% extracellular compartment [e] is merged into the shared pool [u] by
% createMultipleSpeciesModel function. The synaptic cleft [syn] is also merged
% into a shared synaptic pool [q] via reversible transfer reactions, enabling metabolite
% exchange across the synaptic cleft.
%
% Output: iMMU1868-TPS.mat
% ============================================================================

clc; clear;

cd(fileparts(which(mfilename)));

%% Toolbox setup
% Edit toolbox path
% COBRA Toolbox : https://opencobra.github.io/cobratoolbox/
% Gurobi        : https://www.gurobi.com
COBRA_PATH = 'C:/cobra';                  
SOLVER_PATH = 'C:/gurobi1202/win64/matlab';

restoredefaultpath; rehash toolboxcache;
addpath(genpath(COBRA_PATH));
addpath(SOLVER_PATH);

initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');

%% 1) Input files (expanded cell-line specific models)
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end

astro_file = fullfile(scriptDir, 'iMMU1868a-pap.xml');
preNeuron_file = fullfile(scriptDir, 'iMMU1868n-pre.xml');
postNeuron_file = fullfile(scriptDir, 'iMMU1868n-post.xml');

outDir = scriptDir;

%% 2) Load models
fprintf('Loading models...\n');

AST = readCbModel(astro_file);
PRE = readCbModel(preNeuron_file);
POST = readCbModel(postNeuron_file);

% Confirm each model loaded as a valid COBRA struct
assert(isstruct(AST)  && isfield(AST,  'S') && isfield(AST,  'rxns'), ...
       'Astrocyte model failed to load as COBRA struct');
assert(isstruct(PRE)  && isfield(PRE,  'S') && isfield(PRE,  'rxns'), ...
       'Pre-neuron model failed to load as COBRA struct');
assert(isstruct(POST) && isfield(POST, 'S') && isfield(POST, 'rxns'), ...
       'Post-neuron model failed to load as COBRA struct');

% Check the composition of each loaded model
fprintf('  Astrocyte (PAP):   rxns=%d, mets=%d, genes=%d, comps=%d\n', ...
    numel(AST.rxns), numel(AST.mets), numel(AST.genes), numel(AST.comps));
fprintf('  Pre-neuron:        rxns=%d, mets=%d, genes=%d, comps=%d\n', ...
    numel(PRE.rxns), numel(PRE.mets), numel(PRE.genes), numel(PRE.comps));
fprintf('  Post-neuron:       rxns=%d, mets=%d, genes=%d, comps=%d\n', ...
    numel(POST.rxns), numel(POST.mets), numel(POST.genes), numel(POST.comps));

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
fprintf('\nBuilding community model with tags: %s, %s, %s\n', ...
    desiredTags{1}, desiredTags{2}, desiredTags{3});

modelCom = createMultipleSpeciesModel({AST; PRE; POST}, biomassNames, ...
    'nameTagsModels', desiredTags);

modelCom.csense = repmat('E', 1, numel(modelCom.mets));

detectedTags = detectSpeciesPrefixes(modelCom.rxns, 3);
[modelCom.infoCom, modelCom.indCom] = getMultiSpeciesModelId(modelCom, detectedTags);

fprintf('Detected rxn prefixes: %s, %s, %s\n', ...
    detectedTags{1}, detectedTags{2}, detectedTags{3});
fprintf('Community: rxns=%d, mets=%d\n', ...
    numel(modelCom.rxns), numel(modelCom.mets));

%% 5) Remove EX_ "species" if present
if isfield(modelCom.infoCom, 'spAbbr')
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
fprintf('Species count: %d\n', numel(modelCom.infoCom.spAbbr));

%% 6) Rebuild gene fields from individual species models
if ~isfield(modelCom,'genes') || isempty(modelCom.genes) || ...
   ~isfield(modelCom,'rxnGeneMat') || size(modelCom.rxnGeneMat,2)==0
    modelCom = restoreGeneFieldsFromSpecies(modelCom, {AST, PRE, POST}, detectedTags);
    fprintf('Rebuilt gene fields: total genes = %d\n', numel(modelCom.genes));
end

%% 7) Map species exchange bounds into community [u] exchanges
mapToU  = @(metIDs) strrep(metIDs, '[e]', '[u]');
getExchMet = @(mod) find(~cellfun('isempty', regexp(mod.mets, '\[e\]')));

metEx_AST  = getExchMet(AST);  rxn_AST  = mapE_MetToExRxn(AST,  metEx_AST);
metEx_PRE  = getExchMet(PRE);  rxn_PRE  = mapE_MetToExRxn(PRE,  metEx_PRE);
metEx_POST = getExchMet(POST);  rxn_POST = mapE_MetToExRxn(POST, metEx_POST);

lb_AST = nan(numel(metEx_AST),1);
lb_PRE = nan(numel(metEx_PRE),1);
lb_POST = nan(numel(metEx_POST),1);

lb_AST(~isnan(rxn_AST)) = AST.lb(rxn_AST(~isnan(rxn_AST)));
lb_PRE(~isnan(rxn_PRE)) = PRE.lb(rxn_PRE(~isnan(rxn_PRE)));
lb_POST(~isnan(rxn_POST)) = POST.lb(rxn_POST(~isnan(rxn_POST)));

[found_AST, idxComMet_AST] = ismember(mapToU(AST.mets(metEx_AST)),   modelCom.infoCom.Mcom);
[found_PRE, idxComMet_PRE] = ismember(mapToU(PRE.mets(metEx_PRE)),   modelCom.infoCom.Mcom);
[found_POST, idxComMet_POST] = ismember(mapToU(POST.mets(metEx_POST)), modelCom.infoCom.Mcom);

modelCom.lb(modelCom.indCom.EXcom(:,1)) = 0;
modelCom.ub(modelCom.indCom.EXcom(:,1)) = 1000;

for t = find(found_AST).'
    u = idxComMet_AST(t);
    j = modelCom.indCom.EXcom(u,1);
    modelCom.lb(j) = min(modelCom.lb(j), lb_AST(t));
end

for t = find(found_PRE).'
    u = idxComMet_PRE(t);
    j = modelCom.indCom.EXcom(u,1);
    modelCom.lb(j) = min(modelCom.lb(j), lb_PRE(t));
end

for t = find(found_POST).'
    u = idxComMet_POST(t);
    j = modelCom.indCom.EXcom(u,1);
    modelCom.lb(j) = min(modelCom.lb(j), lb_POST(t));
end

fprintf('Community exchange bounds mapped from species models.\n');

%% 7b) Refresh infoCom/indCom before adding the shared pool
[modelCom.infoCom, modelCom.indCom] = getMultiSpeciesModelId(modelCom, detectedTags);
fprintf('infoCom/indCom refreshed (EXcom: %d, EXsp: %d x %d)\n', ...
    size(modelCom.indCom.EXcom,1), size(modelCom.indCom.EXsp,1), size(modelCom.indCom.EXsp,2));

%% 7c) Remove unintended EX_ species
if isfield(modelCom.infoCom, 'spAbbr')
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
        fprintf('Removed spurious EX_ species after infoCom refresh.\n');
    end
end
fprintf('Species count after cleanup: %d\n', numel(modelCom.infoCom.spAbbr));

%% 8) Add the shared synaptic cleft pool [q]
cleftComp = 'syn';   % compartment lable in individual models
cleftPoolComp = 'q';     % shared pool label for the cleft

modelCom = addSharedPoolForCompartment(modelCom, detectedTags, ...
    cleftComp, cleftPoolComp);

fprintf('\n Synaptic Cleft Shared Pool Summary\n');
synPoolMets = modelCom.mets(contains(modelCom.mets, ['[' cleftPoolComp ']']));
fprintf('Shared pool [%s] metabolites: %d\n', cleftPoolComp, numel(synPoolMets));
for k = 1:numel(synPoolMets)
    fprintf('  %s\n', synPoolMets{k});
end

%% 9) Community model summary
fprintf('\n COMMUNITY MODEL SUMMARY: \n');
fprintf('Total reactions: %d\n', numel(modelCom.rxns));
fprintf('Total metabolites: %d\n', numel(modelCom.mets));
fprintf('Total genes: %d\n', numel(modelCom.genes));
fprintf('Species: %s\n', strjoin(modelCom.infoCom.spAbbr, ', '));
fprintf('EXcom (environment exchanges): %d\n', size(modelCom.indCom.EXcom, 1));
fprintf('EXsp  (species exchanges): %d x %d\n', ...
    size(modelCom.indCom.EXsp, 1), size(modelCom.indCom.EXsp, 2));

icleft_rxns = modelCom.rxns(startsWith(modelCom.rxns, 'ICLEFT_'));
fprintf('Synaptic coupling transfers:   %d\n', numel(icleft_rxns));

%% 10) Save community model
matFile = fullfile(outDir, 'iMMU1868-TPS.mat');
if exist(matFile, 'file')
    fileattrib(matFile, '+w', 'a');
    delete(matFile);
end
save(matFile, 'modelCom', '-v7.3');
fprintf('\nSaved community model: %s\n', matFile);

%% Local Functions
% Required by the script to run properly. 

function rxnIdx = mapE_MetToExRxn(model, extracellularMets)
% For each extracellular metabolite, find its exchange reaction (NaN if none).
    singleMetRxns = find(sum(model.S ~= 0, 1) == 1);   % reactions with one metabolite

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
    poolBracket = ['[', compPool, ']'];

    addedPoolMets = containers.Map;

    for si = 1:numel(tags)
        tag = tags{si};

        speciesMetIdx = find(startsWith(modelCom.mets, tag) & ...
                        ~cellfun('isempty', regexp(modelCom.mets, fromPattern)));

        for mi = 1:numel(speciesMetIdx)
            speciesMet = modelCom.mets{speciesMetIdx(mi)};

            metNoTag = erase(speciesMet, tag);
            poolMet  = regexprep(metNoTag, fromPattern, poolBracket);

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
                        'lowerBound', 0, ...
                        'upperBound', 0);
                else
                    modelCom = changeRxnBounds(modelCom, exRxn, 0, 'l');
                    modelCom = changeRxnBounds(modelCom, exRxn, 0, 'u');
                end
            end

            transferRxn = ['ICLEFT_', tag, ...
                     regexprep(poolMet, ['\[', compPool, '\]'], ''), ...
                     '_', compFrom, '_to_', compPool];

            if findRxnIDs(modelCom, transferRxn) == 0
                modelCom = addReaction(modelCom, transferRxn, ...
                    'metaboliteList', {speciesMet, poolMet}, ...
                    'stoichCoeffList', [1, -1], ...
                    'lowerBound', -1000, ...
                    'upperBound', 1000);
            end
        end
    end

    fprintf('Added shared cleft pool [%s] from species compartment [%s]\n', ...
        compPool, compFrom);
end


function modelCom = restoreGeneFieldsFromSpecies(modelCom, speciesModels, tags)
% Rebuild genes, rxnGeneMat, rules, and grRules (with species prefixes)
% from the individual models, since createMultipleSpeciesModel can drop them.
    nRxnsCom    = numel(modelCom.rxns);
    nSpecies    = numel(speciesModels);
    nGenesPerSp = zeros(nSpecies, 1);

    for si = 1:nSpecies
        spModel = speciesModels{si};
        if isfield(spModel, 'genes') && ~isempty(spModel.genes)
            nGenesPerSp(si) = numel(spModel.genes);
        else
            nGenesPerSp(si) = 0;
            spModel.genes      = {};
            spModel.rxnGeneMat = sparse(numel(spModel.rxns), 0);
            if ~isfield(spModel, 'rules'),   spModel.rules   = repmat({''}, numel(spModel.rxns), 1); end
            if ~isfield(spModel, 'grRules'), spModel.grRules = repmat({''}, numel(spModel.rxns), 1); end
            speciesModels{si} = spModel;
        end
    end

    nGenesTotal = sum(nGenesPerSp);
    allGenes = cell(nGenesTotal, 1);
    allGeneNames = cell(nGenesTotal, 1);
    geneColStart = [0; cumsum(nGenesPerSp(1:end-1))] + 1;

    rxnGeneMatCom = sparse(nRxnsCom, nGenesTotal);
    rules   = repmat({''}, nRxnsCom, 1);
    grRules = repmat({''}, nRxnsCom, 1);

    comRxnMap = containers.Map(modelCom.rxns, num2cell(1:nRxnsCom));

    for si = 1:nSpecies
        spModel = speciesModels{si};
        tag  = tags{si};
        gStart = geneColStart(si);
        nGenes = nGenesPerSp(si);

        for gi = 1:nGenes
            allGenes{gStart+gi-1} = [tag, spModel.genes{gi}];
            if isfield(spModel, 'geneNames') && numel(spModel.geneNames) >= gi && ...
               ~isempty(spModel.geneNames{gi})
                allGeneNames{gStart+gi-1} = spModel.geneNames{gi};
            else
                allGeneNames{gStart+gi-1} = spModel.genes{gi};
            end
        end

        % One combined regex to prefix all genes at once.
        if nGenes > 0
            escapedGenes  = cellfun(@regexpescape, spModel.genes, 'UniformOutput', false);
            [~, lenOrder] = sort(cellfun(@numel, spModel.genes), 'descend');
            sortedEscaped = escapedGenes(lenOrder);
            genePattern   = ['(?<![A-Za-z0-9_])(?:', strjoin(sortedEscaped, '|'), ')(?![A-Za-z0-9_])'];
        else
            genePattern = '';
        end

        geneIdxOffset = gStart - 1;

        for ri = 1:numel(spModel.rxns)
            comRxnID = [tag, spModel.rxns{ri}];
            if ~isKey(comRxnMap, comRxnID), continue; end
            comRxnIdx = comRxnMap(comRxnID);

            if nGenes > 0 && isfield(spModel, 'rxnGeneMat') && ~isempty(spModel.rxnGeneMat)
                geneCols = find(spModel.rxnGeneMat(ri,:));
                if ~isempty(geneCols)
                    rxnGeneMatCom(comRxnIdx, gStart-1+geneCols) = spModel.rxnGeneMat(ri, geneCols);
                end
            end

            if isfield(spModel, 'rules') && numel(spModel.rules) >= ri && ~isempty(spModel.rules{ri})
                rules{comRxnIdx} = remapRuleIndices(spModel.rules{ri}, geneIdxOffset);
            end

            if isfield(spModel, 'grRules') && numel(spModel.grRules) >= ri && ~isempty(spModel.grRules{ri})
                if ~isempty(genePattern)
                    grRules{comRxnIdx} = regexprep(spModel.grRules{ri}, genePattern, [tag, '$0']);
                else
                    grRules{comRxnIdx} = spModel.grRules{ri};
                end
            end
        end
        fprintf('  Gene mapping for species %s complete (%d rxns)\n', tag, numel(spModel.rxns));
    end

    modelCom.genes = allGenes;
    modelCom.geneNames = allGeneNames;
    modelCom.rxnGeneMat = rxnGeneMatCom;

    if any(~cellfun(@isempty, rules)),   modelCom.rules   = rules;   end
    if any(~cellfun(@isempty, grRules)), modelCom.grRules = grRules; end
end


function out = remapRuleIndices(ruleIn, offset)
% Shift each gene index x(N) to x(N+offset) for community indexing.
    if offset == 0
        out = ruleIn;
        return;
    end
    idxTokens = regexp(ruleIn, 'x\((\d+)\)', 'tokens');
    tokStart = regexp(ruleIn, 'x\((\d+)\)', 'start');
    tokEnd = regexp(ruleIn, 'x\((\d+)\)', 'end');
    if isempty(idxTokens)
        out = ruleIn;
        return;
    end
    out = '';
    prevEnd = 0;
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
    for ch = metaChars
        str = strrep(str, ch, ['\', ch]);
    end
end
