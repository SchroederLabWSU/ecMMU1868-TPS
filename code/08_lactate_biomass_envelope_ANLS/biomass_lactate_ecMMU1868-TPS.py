"""
Biomass-Lactate Production Envelope for ecMMU1868-TPS (enzyme-constrained community)

This script computes the biomass-lactate production envelope by determining
the minimum and maximum feasible lactate exchange rates across the full range
of community growth rates.

Equal-growth coupling and prot_pool floor (lb = -10) are already baked into
ecMMU1868-TPS.mat by Tripartite_ecCom_equalGrowth.m.

You need to copy the .mat model from Zenodo to run this locally,
since the model file is too large for GitHub. Or run the .m files in
code/07_steadyCom_community/02_ec_TPS_community to generate it locally.

Input:  models/06_community_TPS/enzyme_constrained_tps/ecMMU1868-TPS.mat
Output: results/08_lactate_biomass_envelope_ANLS/biomass_lactate_ecMMU1868-TPS.xlsx
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
RESULTS = REPO_ROOT / "results" / "08_lactate_biomass_envelope_ANLS"
MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
OUTPUT_XLSX = RESULTS / "biomass_lactate_ecMMU1868-TPS.xlsx"

sys.path.insert(0, str(SCRIPT_DIR))
from mat73_to_cobra import load_mat73_cobra

N_PTS = 21 # number of points in the biomass-lactate envelope
TOL = 1e-6 # tolerance used when fixing biomass flux


print(f"Loading {MAT_MODEL} ...", flush=True)
model = load_mat73_cobra(str(MAT_MODEL))
for solver in ("gurobi", "glpk"):
    try: model.solver = solver; break
    except: continue
print(f"  {len(model.reactions)} rxns, {len(model.metabolites)} mets", flush=True)

# Equal-growth coupling, medium, and prot_pool floor are already baked into ecMMU1868-TPS.mat
ast_biomass = model.reactions.get_by_id("AST_BIOMASS_reaction")
lactate_rxn = model.reactions.get_by_id("EX_lac__L[u]")

# Any of the cell biomass reaction can be used because equal-growth coupling enforces v_AST = v_PRE = v_POST.
model.objective = ast_biomass
mu_max = model.optimize().objective_value # maximum community growth rate
print(f"mu_max = {mu_max:.6f} h-1 ", flush=True)
 
rows = []
bm_lb0, bm_ub0 = ast_biomass.lower_bound, ast_biomass.upper_bound

# Fix biomass at fractions of mu_max and determine the feasible lactate secretion range at each growth level.
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
