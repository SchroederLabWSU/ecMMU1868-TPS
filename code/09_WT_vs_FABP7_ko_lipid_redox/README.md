# FABP7 Knockout: Lipid and Redox Analysis

This analysis quantifies FABP7-dependent clearance of neuron-derived peroxidized lipid (13-HPODE; C04717) in the enzyme-constrained tripartite synapse model. Additional analyses test the sensitivity of the clearance phenotype to FABP5/FABP3 turnover numbers, metabolic objective, and alternative optimal flux distributions.

Peroxidized lipid is delivered at the presynaptic neuron. Clearance is measured downstream of FABP transport through astrocytic glutathione-dependent reduction.

Inputs: community model in models/06_community_TPS/enzyme_constrained_tps (see note on large files in the main README).

Files:
1. 1_kcat_load_sweep.py
Computes wild-type and FABP7-knockout clearance across increasing delivered 13-HPODE loads. FABP5 and FABP3 turnover numbers are varied from 0.1x to 10x their DLKcat predictions while the wild-type reference is held at the baseline turnover values.
2. 2_objective_comparison.py
Tests whether the clearance phenotype changes with ATPM demand. Clearance is evaluated under the biomass baseline, with biomass blocked and ATPM free, and at increasing fractions of maximal ATPM demand. The biomass baseline and 0.99 ATPM condition are saved, while the remaining conditions are printed as a robustness check.
3. 3_fva_clearance.py
Computes the minimum and maximum feasible clearance flux in wild type and FABP7 knockout while the objective remains optimal. This tests whether the clearance phenotype persists across alternative optimal flux distributions.
4. mat73_to_cobra.py
Is the helper function for loading MATLAB v7.3 community models into COBRApy.


Outputs: Results go to results/09_WT_vs_FABP7_ko_lipid_redox:
1. kcat_load_sweep.xlsx - FABP parameters, clearance across delivered loads, and FABP5/FABP3 turnover sensitivity.
2. objective_comparison.xlsx - clearance under the biomass and ATPM conditions.
3. fva_clearance.xlsx - minimum and maximum feasible clearance flux for wild type and FABP7 knockout.
