"""
Biomass-Lactate Production Envelope for iMMU1868-TPS (stoichiometric TPS community)

Computes the biomass-lactate production envelope by determining the minimum
and maximum feasible lactate secretion rates across the full range of
community growth rates. This is the stoichiometric (no enzyme constraints) version of the analysis.

Equal-growth coupling (v_AST = v_PRE = v_POST) is added here because it is
not already included in iMMU1868-TPS.mat.

You need to copy the .mat model from Zenodo to run this locally,
since the model file is too large for GitHub. Or run Tripartite_Model_SteadyCom.m
in code/07_steadyCom_community/01_stoichiometric_TPS_community to generate it locally.

Input:  models/06_community_TPS/stoichiometric_tps/iMMU1868-TPS.mat
Output: results/08_lactate_biomass_envelope_ANLS/biomass_lactate_iMMU1868-TPS.xlsx
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
RESULTS = REPO_ROOT / "results" / "08_lactate_biomass_envelope_ANLS"
MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "stoichiometric_tps" / "iMMU1868-TPS.mat"
OUTPUT_XLSX = RESULTS / "biomass_lactate_iMMU1868-TPS.xlsx"

sys.path.insert(0, str(SCRIPT_DIR))
from mat73_to_cobra import load_mat73_cobra

N_PTS = 21 # number of points in the biomass-lactate envelope
TOL = 1e-6 # tolerance used when fixing biomass flux

# Defined minimal medium for the community [u] pool
MEDIUM = [
    ("EX_h2o[u]", -1000, 1000), ("EX_h[u]", -1000, 1000),
    ("EX_na1[u]", -1000, 1000), ("EX_k[u]", -1000, 1000),
    ("EX_cl[u]", -1000, 1000), ("EX_ca2[u]", -1000, 1000),
    ("EX_nh4[u]", -10, 1000), ("EX_pi[u]", -10, 1000),
    ("EX_so4[u]", -10, 1000), ("EX_o2[u]", -10, 1000),
    ("EX_co2[u]",  0, 1000), ("EX_glc__D[u]", -10, 1000),
    ("EX_gln__L[u]", 0, 1000), ("EX_lac__L[u]", 0, 1000),
    ("EX_hdca[u]", 0, 1000), ("EX_val__L[u]", -0.155172, 1000),
    ("EX_his__L[u]", -0.186207, 1000), ("EX_ile__L[u]", -0.186207, 1000),
    ("EX_leu__L[u]", -0.186207, 1000), ("EX_lys__L[u]", -0.186207, 1000),
    ("EX_met__L[u]", -0.155172, 1000), ("EX_phe__L[u]", -0.279310, 1000),
    ("EX_thr__L[u]", -0.124138, 1000), ("EX_trp__L[u]", -0.341379, 1000),
]


def add_equal_growth_coupling(model):
    # Force v_AST = v_PRE = v_POST via two equality constraints.
    ast_biomass = model.reactions.get_by_id("AST_BIOMASS_reaction")
    pre_biomass = model.reactions.get_by_id("PRE_BIOMASS_reaction")
    post_biomass = model.reactions.get_by_id("POST_BIOMASS_reaction")
    model.add_cons_vars([
        model.problem.Constraint(ast_biomass.flux_expression - pre_biomass.flux_expression,
              lb=0, ub=0, name="couple_AST_PRE"),
        model.problem.Constraint(ast_biomass.flux_expression - post_biomass.flux_expression,
              lb=0, ub=0, name="couple_AST_POST"),
    ])


def apply_medium(model, medium):
    #Close all community exchanges then open the defined medium.
    for rxn in model.reactions:
        if rxn.id.startswith("EX_") and rxn.id.endswith("[u]"):
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
    rxn_set = {r.id for r in model.reactions}
    for rxn_id, lb, ub in medium:
        if rxn_id in rxn_set:
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.lower_bound = lb
            rxn.upper_bound = ub


print(f"Loading {MAT_MODEL} .....", flush=True)
model = load_mat73_cobra(str(MAT_MODEL))
for solver in ("gurobi", "glpk"):
    try: model.solver = solver; break
    except: continue
print(f"  {len(model.reactions)} rxns, {len(model.metabolites)} mets", flush=True)

add_equal_growth_coupling(model)
apply_medium(model, MEDIUM)

ast_biomass = model.reactions.get_by_id("AST_BIOMASS_reaction")
lactate_rxn = model.reactions.get_by_id("EX_lac__L[u]")

# Any biomass reaction can be used because equal-growth coupling enforces v_AST = v_PRE = v_POST.
model.objective = ast_biomass
mu_max = model.optimize().objective_value
print(f"mu_max = {mu_max:.6f} h-1", flush=True)

rows = []
bm_lb0, bm_ub0 = ast_biomass.lower_bound, ast_biomass.upper_bound
for i, fraction in enumerate(np.linspace(0.0, 1.0, N_PTS), 1):
    biomass = fraction * mu_max
    ast_biomass.lower_bound = biomass - TOL
    ast_biomass.upper_bound = biomass + TOL

    model.objective = lactate_rxn
    model.objective_direction = "max"
    sol_max = model.optimize()
    lac_max = sol_max.objective_value if sol_max.status == "optimal" else float("nan")

    model.objective_direction = "min"
    sol_min = model.optimize()
    lac_min = sol_min.objective_value if sol_min.status == "optimal" else float("nan")
    
    # Restore default biomass bounds before the next envelope point.
    ast_biomass.lower_bound = bm_lb0
    ast_biomass.upper_bound = bm_ub0
    rows.append({"fraction_muMax": fraction, "biomass_flux_h-1": biomass,
                 "lac_min_mmol_gDW_h": lac_min, "lac_max_mmol_gDW_h": lac_max})
    print(f"  [{i:2d}/{N_PTS}] frac={fraction:.2f}  bm={biomass:.4f}  "
          f"lac_min={lac_min:.4f}  lac_max={lac_max:.4f}", flush=True)

df = pd.DataFrame(rows).sort_values("biomass_flux_h-1").reset_index(drop=True)

# Save results
RESULTS.mkdir(parents=True, exist_ok=True)
df.to_excel(OUTPUT_XLSX, index=False)
print(f"\nSaved: {OUTPUT_XLSX}")
