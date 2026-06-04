"""
C04717 Input Sensitivity Analysis on ecMMU1868-TPS

This determines the maximum feasible peroxidized-linoleate (13-HPODE, C04717)
input that can be sustained by the enzyme-constrained tripartite synapse
community model under WT and FABP7-KO conditions.

Output: results/09_WT_vs_FABP7_ko_lipid_redox4/C04717_input_sweep.xlsx
"""

import sys
from pathlib import Path
import pandas as pd
from cobra import Reaction

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
RESULTS = REPO_ROOT / "results" / "09_WT_vs_FABP7_ko_lipid_redox4"
OUTPUT_XLSX = RESULTS / "C04717_input_sweep.xlsx"

sys.path.insert(0, str(SCRIPT_DIR))
from mat73_to_cobra import load_mat73_cobra

MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"

FABP7_GENE = "12140" # Entrez gene ID for FABP7
EPS_DETOX = 1e-4 # secondary objective weight on astrocyte lipid-repair reaction (AST_LNLCOH_REPAIRc)

# Explicit reaction IDs in the community model
HPODE_SOURCE = "PRE_SK_C04717_prec"
HPODE_PREC_METS = ["PRE_C04717[prec]", "PRE_C04717_prec"]
HPODE_TRANSPORT = "PRE_C04717tprec"
LINOLEATE_SINK = "AST_SK_lnlc_c"
LIPID_REPAIR = "AST_LNLCOH_REPAIRc"
BIOMASS_AST = "AST_BIOMASS_reaction"

FABP7_KO_EXTRA = ["AST_C04717tc", "AST_C04717td", "AST_C04717td_REV", "AST_C04717ATP"]
INPUT_SWEEP_GRID = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0]


def ensure_src(model):
    # Add the 13-HPODE source reaction if not already in the model.
    if HPODE_SOURCE in {r.id for r in model.reactions}:
        return
    hpode_met = next((c for c in HPODE_PREC_METS if c in {m.id for m in model.metabolites}), None)
    assert hpode_met, "PRE C04717 prec metabolite not found"
    src = Reaction(HPODE_SOURCE)
    src.bounds = (0, 0)
    src.add_metabolites({model.metabolites.get_by_id(hpode_met): 1.0})
    model.add_reactions([src])


def set_force(model, flux_value):
    # Set the forced 13-HPODE input to flux_value.
    model.reactions.get_by_id(HPODE_SOURCE).bounds = (flux_value, flux_value)
    transport = model.reactions.get_by_id(HPODE_TRANSPORT)
    transport.bounds = (max(0.0, flux_value), max(transport.upper_bound, flux_value, 1000.0))


def solve_at(model, flux_value, label):
    # Solve at a given C04717 input and return (feasible, biomass, repair).
    set_force(model, flux_value)
    sol = model.optimize()
    feasible   = sol.status == "optimal"
    biomass    = sol.fluxes.get(BIOMASS_AST) if feasible else float("nan")
    repair_flux= sol.fluxes.get(LIPID_REPAIR) if feasible else float("nan")
    print(f"  [{label}] C04717={flux_value:>8.3g}  {sol.status:<12}  "
          f"AST_BM={biomass if feasible else 'NA'}  repair={repair_flux if feasible else 'NA'}", flush=True)
    return feasible, biomass, repair_flux


def binary_search_max(model, label, lo, hi, tol=0.05):
    # Find the largest feasible C04717 input via binary search.
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        feasible, _, _ = solve_at(model, mid, f"{label} bsearch")
        lo, hi = (mid, hi) if not feasible else (lo, mid)
    return lo if not solve_at(model, hi, label)[0] else hi


print(f"Loading {MAT_MODEL} .....", flush=True)
m = load_mat73_cobra(str(MAT_MODEL))
for solver in ("gurobi", "glpk"):
    try: m.solver = solver; break
    except: continue
print(f"  {len(m.reactions)} rxns", flush=True)

ensure_src(m)
m.reactions.get_by_id(LINOLEATE_SINK).upper_bound = 1000

ast_biomass  = m.reactions.get_by_id(BIOMASS_AST)
lipid_repair = m.reactions.get_by_id(LIPID_REPAIR)
m.objective  = {ast_biomass: 1.0, lipid_repair: EPS_DETOX}

rows = []

# WT sweep 
print("\nWT coarse sweep:", flush=True)
wt_last_feasible = None; wt_first_infeasible = None
for x in INPUT_SWEEP_GRID:
    feasible, biomass, repair_flux = solve_at(m, x, "WT")
    rows.append({"condition": "WT", "C04717_input": x,
         "status": "optimal" if feasible else "infeasible",
         "AST_BIOMASS": biomass, "AST_LNLCOH_REPAIRc": repair_flux})
    if feasible: wt_last_feasible = x
    else: wt_first_infeasible = x; break

if wt_last_feasible and wt_first_infeasible:
    wt_max = binary_search_max(m, "WT", wt_last_feasible, wt_first_infeasible)
    feasible, biomass, repair_flux = solve_at(m, wt_max, "WT_MAX")
    rows.append({"condition": "WT_MAX", "C04717_input": wt_max,
            "status": "optimal" if feasible else "infeasible",
            "AST_BIOMASS": biomass, "AST_LNLCOH_REPAIRc": repair_flux})
    print(f"WT max feasible C04717 = {wt_max:.4g} mmol/gDW/h", flush=True)

# FABP7 KO sweep
print("\nFABP7 KO sweep:", flush=True)
ko_blocked = [r.id for r in m.reactions
         if r.gene_reaction_rule and FABP7_GENE in r.gene_reaction_rule]
for rxn_id in FABP7_KO_EXTRA:
    if rxn_id in {rr.id for rr in m.reactions} and rxn_id not in ko_blocked:
        ko_blocked.append(rxn_id)
for rxn_id in ko_blocked:
    m.reactions.get_by_id(rxn_id).bounds = (0, 0)
print(f"  Closed {len(ko_blocked)} reactions for FABP7 KO", flush=True)

ko_last_feasible = None; ko_first_infeasible = None
for x in INPUT_SWEEP_GRID:
    feasible, biomass, repair_flux = solve_at(m, x, "KO")
    rows.append({"condition": "KO", "C04717_input": x,
            "status": "optimal" if feasible else "infeasible",
            "AST_BIOMASS": biomass, "AST_LNLCOH_REPAIRc": repair_flux})
    if feasible: ko_last_feasible = x
    else: ko_first_infeasible = x; break

if ko_last_feasible and ko_first_infeasible:
    ko_max = binary_search_max(m, "KO", ko_last_feasible, ko_first_infeasible)
    feasible, biomass, repair_flux = solve_at(m, ko_max, "KO_MAX")
    rows.append({"condition": "KO_MAX", "C04717_input": ko_max,
            "status": "optimal" if feasible else "infeasible",
            "AST_BIOMASS": biomass, "AST_LNLCOH_REPAIRc": repair_flux})
    print(f"KO max feasible C04717 = {ko_max:.4g} mmol/gDW/h", flush=True)

# Save results
RESULTS.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_excel(OUTPUT_XLSX, index=False)
print(f"\nSaved: {OUTPUT_XLSX}")
