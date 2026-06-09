# Enzyme-Constrained TPS Community (ecMMU1868-TPS)

This builds the enzyme-constrained tripartite synapse community model, ecMMU1868-TPS, from the three enzyme-constrained cell models:

1. ecMMU1868a-pap
2. ecMMU1868n-pre
3. ecMMU1868n-post

The models are combined using a modified SteadyCom framework. This differs from standard SteadyCom because each cell model contains GECKO protein-pool pseudoreactions, and the TPS model includes a shared synaptic-cleft exchange pool plus the standard shared extracellular pool.

Inputs: from the folder models/05_enzyme_constrained

Outputs: Saved in models/06_community_TPS/enzyme_constrained_tps

1. ecTPS_build.mat
   The enzyme-constrained community model before constraints are applied.

2. ecMMU1868-TPS.mat
   Final enzyme-constrained TPS community model with equal-growth coupling, defined medium, protein-pool bound, and FBA objective baked in. This is the model used in all downstream analyses. ecTPS_build.mat can also be used for analysis, but the equal-growth coupling must be applied before solving.

Both .mat files exceed GitHub's 100 MB per-file limit and are therefore not stored in this repository. They are archived on Zenodo (https://doi.org/10.5281/zenodo.20606464). To reproduce them, run the two scripts in order.

The Matlab scripts run in order:

1. Tripartite_ecCom.m
   Combines the three ecModels into one SteadyCom community and writes ecTPS_build.mat.

2. Tripartite_ecCom_equalGrowth.m
   Loads ecTPS_build.mat, adds equal-growth coupling, defined medium, protein-pool bound, and FBA objective, then writes ecMMU1868-TPS.mat. The model is saved before the verify step so all annotation fields are preserved.


The final model, ecMMU1868-TPS is used in subsequent analyses:

1. Step 08: biomass-lactate envelope
2. Step 09: FABP7 knockout lipid/redox analysis
3. Step 10: reporter metabolite analysis