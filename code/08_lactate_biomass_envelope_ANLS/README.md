# Biomass-Lactate Envelope (ANLS)

This step evaluates the astrocyte-neuron lactate shuttle (ANLS) by computing biomass-lactate production envelopes for the stoichiometric (iMMU1868-TPS) and enzyme-constrained (ecMMU1868-TPS) tripartite synapse community models. At each growth level, the minimum and maximum feasible lactate secretion rates are determined.

Inputs: community models in models/06_community_TPS (see note on large files below).

Files:
1. biomass_lactate_iMMU1868-TPS.py
Computes the biomass-lactate envelope for the stoichiometric community model. Equal-growth coupling is added within the script because it is not stored in iMMU1868-TPS.mat.
2. biomass_lactate_ecMMU1868-TPS.py
Computes the biomass-lactate envelope for the enzyme-constrained community model. Equal-growth coupling, defined medium, and protein-pool bounds are already incorporated into ecMMU1868-TPS.mat.
3. mat73_to_cobra.py
Is the helper function for loading MATLAB v7.3 community models into COBRApy.


Outputs: Results go to results/08_lactate_biomass_envelope_ANLS:
1. biomass_lactate_iMMU1868-TPS.xlsx - stoichiometric envelope data.
2. biomass_lactate_ecMMU1868-TPS.xlsx - enzyme-constrained envelope data.

Note that the community .mat models exceed the GitHub 100 MB file-size limit and are not stored in this repository. Download them from Zenodo or regenerate them using the scripts in 07_steadyCom_community folder.
