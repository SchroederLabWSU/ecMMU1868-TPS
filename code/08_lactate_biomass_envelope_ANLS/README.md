# Lactate Envelopes (ANLS)

This analysis evaluates the astrocyte-neuron lactate shuttle (ANLS) by computing lactate flux envelopes for the stoichiometric (iMMU1868-TPS) and enzyme-constrained (ecMMU1868-TPS) tripartite synapse community models. At each objective level, the minimum and maximum feasible lactate exchange rates are determined.

Envelopes are calculated under two objectives: community biomass and ATP maintenance with biomass production blocked. ATP maintenance is used as the primary physiological objective, while biomass is retained as a secondary comparison. Additional analyses test the sensitivity of lactate coupling to turnover-number assignment and the protein-pool bound.

Inputs: community models in models/06_community_TPS.

Files:

01_envelopes
1. lactate_envelopes_iMMU1868-TPS.py
Computes biomass- and ATPM-lactate envelopes for the stoichiometric community model.
2. lactate_envelopes_ecMMU1868-TPS.py
Computes biomass- and ATPM-lactate envelopes for the enzyme-constrained community model.

02_kcat_bootstrap
3. kcat_bootstrap.py
Resamples reactions with both BRENDA and DLKcat turnover numbers across 100 parameterizations and recomputes the lactate-coupling onset. Biomass and ATPM analyses use the same random seeds, allowing paired comparison between objectives. Run biomass first, then atpm; the ATPM run reads the biomass workbook to build the paired comparison.

03_protein_pool_sweep
4. protein_pool_sweep.py
Recomputes the biomass-based lactate-coupling onset across 25 protein-pool bounds from 0.1 to 22 mmol/gDW/h.
5. protein_allocation_sweep.py
Quantifies enzyme allocation and protein-budget redistribution across the protein-pool regimes.

6. mat73_to_cobra.py
Is the helper function for loading MATLAB v7.3 community models into COBRApy.


Outputs: Results go to results/08_lactate_envelopes_ANLS:

01_envelopes
1. lactate_envelopes_iMMU1868-TPS.xlsx - stoichiometric envelope data.
2. lactate_envelopes_ecMMU1868-TPS.xlsx - enzyme-constrained envelope data.

02_kcat_bootstrap
3. kcat_bootstrap_biomass.xlsx and kcat_bootstrap_atpm.xlsx - coupling onset across the 100 parameterizations, one file per objective.
4. kcat_parameter_sheet.xlsx - the turnover numbers drawn in each parameterization.

03_protein_pool_sweep
5. protein_pool_sweep.xlsx - lactate-coupling onset across the protein-pool bounds.
6. protein_allocation_sweep.xlsx - enzyme allocation and protein-budget changes across the protein-pool regimes.