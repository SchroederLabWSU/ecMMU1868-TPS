# Model Curation

Curates the published iMM1865 base model into iMMU1867 by fixing stoichiometric inconsistencies, resolving mass and charge imbalances, correcting non-physiological exchange reactions, and standardizing gene-reaction associations and metabolite annotations.

Inputs: models/01_starting_model_iMM1865/iMM1865.xml and models/02_curated_iMMU1867/iMMU1867.xml

Files:
1. model_consistency_check_protocol.ipynb: MEMOTE-style consistency checks on a model.
2. biomassTest_under_defined_medium.ipynb: Confirms the curated model (iMMU1867) grows on the defined minimal medium.
3. iMM1865_vs_iMMU1867_consistency_checks.py: Compares iMM1865 vs iMMU1867 on stoichiometric consistency, mass imbalances, and charge imbalances. Requires memote. 

Outputs: iMM1865_vs_iMMU1867_consistency_checks.py saves results in results/01_model_curation/.