# Enzyme-Constrained TPS Community (ecMMU1868-TPS)

This builds the enzyme-constrained tripartite synapse community model, ecMMU1868-TPS, from the three enzyme-constrained cell models:

1. ecMMU1868a-pap
2. ecMMU1868n-pre
3. ecMMU1868n-post

The models are combined using a modified SteadyCom framework. This differs from standard SteadyCom because each cell model contains GECKO protein-pool pseudoreactions, and the TPS model includes a shared synaptic-cleft exchange pool plus the standard shared extracellular pool.

Inputs: from the folder models/05_enzyme_constrained

Outputs: Saved in models/06_community_TPS/enzyme_constrained_tps

1. ecMMU1868-TPS.mat
   Final enzyme-constrained TPS community model with equal-growth coupling, defined medium, protein-pool bound, and FBA objective baked in. This is the model used in all downstream analyses.

The .mat file exceeds GitHub's 100 MB per-file limit and is therefore not stored in this repository. It is archived on Zenodo (https://doi.org/10.5281/zenodo.20606464). To reproduce it, run the script below.

The Matlab script:

1. Tripartite_ecCom.m
   Combines the three ecModels into one SteadyCom community, then adds equal-growth coupling, defined medium, protein-pool bound, and FBA objective, and writes ecMMU1868-TPS.mat. Community assembly and constraint application are done in a single pass, so no intermediate build file is written.

Note: createMultipleSpeciesModel drops the model.ec metadata field during assembly. The enzyme constraints survive in the S-matrix as 1/kcat coefficients, which is what the downstream analyses use.

The final model, ecMMU1868-TPS is used in subsequent analyses:

1. Step 08: biomass-lactate envelope
2. Step 09: FABP7 knockout lipid/redox analysis
3. Step 10: reporter metabolite analysis
