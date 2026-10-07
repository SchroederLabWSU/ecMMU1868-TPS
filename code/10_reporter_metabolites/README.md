# Reporter Metabolite Analysis

This analysis identifies lipid-associated metabolites whose surrounding reaction fluxes change most between wild type and FABP7 knockout using the Patil-Nielsen reporter-metabolite method. The comparison is repeated across 100 turnover-number parameterizations and evaluated under two objectives: ATP maintenance with growth blocked and community growth.

Reactions with both a BRENDA measurement and a DLKcat prediction are randomly assigned one of the two turnover numbers in each parameterization. The same turnover-number draw is used for both objectives, allowing the robustness of reporter-metabolite rankings to be assessed across parameterizations.

For each reaction, the absolute WT-KO flux difference, |WT - KO|, is standardized across reactions. For each metabolite, the mean standardized change of its connected reactions is scaled by the square root of the number of connected reactions and then standardized within that parameterization. The signed quantity, |WT| - |KO|, is summed across connected reactions to classify the metabolite as suppressed or induced in the knockout.

The full 100-parameterization analysis was run on the Kamiak HPC cluster using SLURM because it requires 400 pFBA solves of the enzyme-constrained community model.

Inputs: community model in models/06_community_TPS/enzyme_constrained_tps (see note on large files in the main README).

Files:
1. 1_extract_metadata.py
Writes the metabolite names, reaction names and S-matrix triplets the scoring needs, so the analysis does not reload the community model.
2. 2_reporter_metabolites.py
Applies each parameterization's turnover-number draw, solves the wild-type and FABP7-knockout models under both objectives, calculates reporter-metabolite scores, ranks lipid-associated metabolites, and writes summary and per-parameterization results. 
3. mat73_to_cobra.py
Is the helper function for loading MATLAB v7.3 community models into COBRApy.

Outputs: Results go to results/10_reporter_metabolites:

1. model_met_meta.csv, model_rxn_meta.csv, model_S_triplets.csv - metabolite names, reaction names, and the S matrix as triplets.
2. kcat_parameters.xlsx - the turnover numbers each reaction can take.
3. reporter_metabolites_bootstrap.xlsx - summary and per-parameterization scores for both objectives.

