"""
WT vs FABP7 Knockout lipid peroxidation simulation in ecMMU1868-TPS.

This script compares WT and FABP7-KO flux distributions in the enzyme-constrained
tripartite synapse community model under forced peroxidized-linoleate input
(C04717, 13-HPODE) into the presynaptic compartment.

The objective maximizes community biomass, with a small secondary objective favoring
astrocytic lipid repair through AST_LNLCOH_REPAIRc.

Outputs in results/09_WT_vs_FABP7_ko_lipid_redox4/:
  full_flux_wt_lipidperox.csv  — full flux distribution for WT
  full_flux_ko_lipidperox.csv  — full flux distribution for FABP7-KO
"""

import sys
from pathlib import Path
import pandas as pd
from cobra import Reaction

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
RESULTS = REPO_ROOT / "results" / "09_WT_vs_FABP7_ko_lipid_redox4"

sys.path.insert(0, str(SCRIPT_DIR))
from mat73_to_cobra import load_mat73_cobra

MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
WT_CSV = RESULTS / "full_flux_wt_lipidperox.csv"
KO_CSV = RESULTS / "full_flux_ko_lipidperox.csv"

FORCE_FLUX = 1e-3  # mmol/gDW/h — forced peroxidized-linoleate (C04717) input
EPS_DETOX = 1e-4  # secondary objective weight on astrocyte lipid-repair reaction (AST_LNLCOH_REPAIRc)
FABP7_GENE = "12140"  # Entrez gene ID for FABP7

# Explicit reaction IDs in the community model
HPODE_SOURCE = "PRE_SK_C04717_prec"  # forces 13-HPODE into pre-synaptic compartment
HPODE_TRANSPORT = "PRE_C04717tprec"  # transports 13-HPODE within pre-synaptic neuron
HPODE_SINK = "PRE_SK_C04717_syn"  # synaptic cleft sink (closed to route HPODE through PAP)
LINOLEATE_SINK = "AST_SK_lnlc_c"  # astrocyte linoleate sink (opened to allow repair)
LIPID_REPAIR = "AST_LNLCOH_REPAIRc"  # astrocyte HODE -> linoleate repair (detox endpoint)
BIOMASS_AST = "AST_BIOMASS_reaction"
BIOMASS_PRE = "PRE_BIOMASS_reaction"
BIOMASS_POST = "POST_BIOMASS_reaction"

# Reactions blocked in FABP7 KO: gene-rule hits + FABP7-dependent bypass routes
FABP7_KO_EXTRA = ["AST_C04717tc", "AST_C04717td", "AST_C04717td_REV", "AST_C04717ATP"]


print(f"Loading {MAT_MODEL} ...", flush=True)
m = load_mat73_cobra(str(MAT_MODEL))
for solver in ("gurobi", "glpk"):
    try: m.solver = solver; break
    except: continue
print(f"  {len(m.reactions)} rxns, {len(m.metabolites)} mets", flush=True)

rxn_set = {r.id for r in m.reactions}

# Force 13-HPODE production into the pre-synaptic compartment
if HPODE_SOURCE in rxn_set:
    m.reactions.get_by_id(HPODE_SOURCE).bounds = (FORCE_FLUX, FORCE_FLUX)
else:
    # Add source reaction if not present in the model
    hpode_met = next((mid for mid in ["PRE_C04717[prec]", "PRE_C04717_prec"]
          if mid in {mm.id for mm in m.metabolites}), None)
    assert hpode_met, "PRE C04717 prec metabolite not found in model"
    src = Reaction(HPODE_SOURCE)
    src.bounds = (FORCE_FLUX, FORCE_FLUX)
    src.add_metabolites({m.metabolites.get_by_id(hpode_met): 1.0})
    m.add_reactions([src])
    print(f"  Added {HPODE_SOURCE} (not in model)", flush=True)

# Ensure transport can carry at least FORCE_FLUX
if HPODE_TRANSPORT in rxn_set:
    tr = m.reactions.get_by_id(HPODE_TRANSPORT)
    tr.lower_bound = max(tr.lower_bound, FORCE_FLUX)

# Close synaptic cleft sink to route HPODE through the astrocyte PAP pathway
if HPODE_SINK in rxn_set:
    m.reactions.get_by_id(HPODE_SINK).bounds = (0, 0)

# Open astrocyte linoleate sink to allow the repair pathway to clear HODE
if LINOLEATE_SINK in rxn_set:
    m.reactions.get_by_id(LINOLEATE_SINK).upper_bound = 1000

# Objective: maximize community biomass with repair as tie-breaker
ast_biomass = m.reactions.get_by_id(BIOMASS_AST)
objective = {ast_biomass: 1.0}
if LIPID_REPAIR in rxn_set:
    objective[m.reactions.get_by_id(LIPID_REPAIR)] = EPS_DETOX
m.objective = objective

# WT simulation: optimize and save fluxes for all reactions in the model
print("\nRunning WT .....", flush=True)
wt_solution = m.optimize()
print(f"  status={wt_solution.status}  AST_BM={wt_solution.fluxes.get(BIOMASS_AST):.6g}", flush=True)
if LIPID_REPAIR in rxn_set:
    print(f"  {LIPID_REPAIR} = {wt_solution.fluxes.get(LIPID_REPAIR):.6g}", flush=True)

all_rxn_ids = [r.id for r in m.reactions]
wt_df = pd.DataFrame({"rxn_id": all_rxn_ids,
                       "flux":   [wt_solution.fluxes.get(rid) for rid in all_rxn_ids]})

# FABP7 KO simulation: block FABP7-associated reactions and re-optimize
print("\nRunning FABP7 KO .....", flush=True)
# Block reactions directly associated with FABP7 and additional FABP7-dependent
# lipid-transfer routes that do not carry explicit gene.
ko_blocked = [r.id for r in m.reactions
    if r.gene_reaction_rule and FABP7_GENE in r.gene_reaction_rule]
for rxn_id in FABP7_KO_EXTRA:
    if rxn_id in rxn_set and rxn_id not in ko_blocked:
        ko_blocked.append(rxn_id)
for rxn_id in ko_blocked:
    m.reactions.get_by_id(rxn_id).bounds = (0, 0)
print(f"  Closed {len(ko_blocked)} reactions for FABP7 KO", flush=True)

ko_solution = m.optimize()
print(f"  status={ko_solution.status}", flush=True)
ko_fluxes = [ko_solution.fluxes.get(rid) if ko_solution.status == "optimal" else float("nan")
     for rid in all_rxn_ids]
ko_df = pd.DataFrame({"rxn_id": all_rxn_ids, "flux": ko_fluxes})

# Save results
RESULTS.mkdir(parents=True, exist_ok=True)
wt_df.to_csv(WT_CSV, index=False)
ko_df.to_csv(KO_CSV, index=False)
print(f"\nSaved: {WT_CSV}")
print(f"Saved: {KO_CSV}")
