classdef ecMMU1868aPapAdapter < ModelAdapter
    % GECKO 3.0 model adapter for the astrocyte+PAP (iMMU1868a-pap) model
    % Compartments: c, e, g, i, l, m, n, pap, r, syn, x.
    
    % Built on the curated mouse model iMMU1867.
    % Cell-type differentiation via iMAT then expanded with the perisynaptic
    % astrocyte process (PAP) and synaptic cleft (syn) compartments
    % for the tripartite synapse.

    methods
        function obj = ecMMU1868aPapAdapter()
            % Path to this model folder (auto-detected via mfilename)
            obj.params.path = fileparts(mfilename('fullpath'));

            % Conventional GEM that this ecModel is built from.
            obj.params.convGEM = fullfile(obj.params.path, 'models', ...
                'iMMU1868a-pap.xml');

            % Average enzyme saturation factor (GECKO default).
            obj.params.sigma = 0.5;

            % Total protein content [g protein / gDW]. GECKO default initial value
            obj.params.Ptot = 0.5;

            % Enzyme mass fraction [g enzyme / g protein]. GECKO default
            % initial value; recalculated from proteomics during the build
            % via calculateFfactor (f = 0.0278 for astrocytes).
            obj.params.f = 0.5;

            % Maximum growth rate the model should be able to reach [1/h]
            obj.params.gR_exp = 0.3269;

            % Organism scientific name
            obj.params.org_name = 'mus musculus';

            % Taxonomic identifier for Complex Portal
            obj.params.complex.taxonomicID = 10090;

            % Organism KEGG ID, selected at
			% https://www.genome.jp/kegg/catalog/org_list.html
            obj.params.kegg.ID = 'mmu';
            obj.params.kegg.geneID = 'kegg';

            % Identifier that should be used to query UniProt.
            % Select proteome IDs at https://www.uniprot.org/proteomes/
            % or taxonomy IDs at https://www.uniprot.org/taxonomy.
            obj.params.uniprot.type = 'taxonomy'; % 'proteome' or 'taxonomy'
            obj.params.uniprot.ID = '10090'; % should match the ID type

            % Field for Uniprot gene ID - should match the gene ids used in the 
            % model. It should be one of the "Returned Field" entries under
            % "Names & Taxonomy" at this page: https://www.uniprot.org/help/return_fields
            obj.params.uniprot.geneIDfield = 'gene_names';

            % Whether only reviewed data from UniProt should be considered.
            % Reviewed data has highest confidence, but coverage might be (very)
            % low for non-model organisms
            obj.params.uniprot.reviewed = false;

            % Reaction ID for glucose exchange reaction
            obj.params.c_source = 'EX_glc__D_e';

            % Reaction ID for biomass pseudoreaction
            obj.params.bioRxn   = 'BIOMASS_reaction';

            % Name of the compartment where the protein pseudometabolites
            % should be located (all be located in the same compartment,
            % this does not interfere with them catalyzing reactions in
            % different compartments). Typically, cytoplasm is chosen.
            obj.params.enzyme_comp = 'Cytoplasm';
        end

        function [spont, spontRxnNames] = getSpontaneousReactions(obj, model)
            % Indicates how spontaneous reactions are identified. Here it
            % is done by the reaction have 'spontaneous' in its name.
            spont = contains(model.rxnNames, 'spontaneous');
            spontRxnNames = model.rxnNames(spont);
        end
    end
end
