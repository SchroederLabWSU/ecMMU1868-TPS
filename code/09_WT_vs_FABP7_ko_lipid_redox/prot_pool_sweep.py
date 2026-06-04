"""
Sensitivity analysis of the protein-pool exchange lower bound in the
enzyme-constrained tripartite synapse community model.

The community growth is evaluated across prot_pool_exchange lower bounds
from -1000 to -0.001 mmol/gDW/h. This was used to select the
final bound of -10 mmol/gDW/h.

The sweep uses ecTPS_build.mat (the reconstructed community model before
equal-growth coupling and final constraints are applied), allowing the
protein-pool bound to be varied freely before the construction of the final
ecMMU1868-TPS model.

Input:  ecTPS_build.mat  (in models/06_community_TPS/enzyme_constrained_tps/)
Output: results/09_WT_vs_FABP7_ko_lipid_redox4/prot_pool_sweep.xlsx
"""

import sys
from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
RESULTS = REPO_ROOT / "results" / "09_WT_vs_FABP7_ko_lipid_redox4"
OUTPUT_XLSX = RESULTS / "prot_pool_sweep.xlsx"

sys.path.insert(0, str(SCRIPT_DIR))
from mat73_to_cobra import load_mat73_cobra

MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecTPS_build.mat"

BIOMASS_ID = "AST_BIOMASS_reaction"
PROT_RXNS = ["AST_prot_pool_exchange", "PRE_prot_pool_exchange", "POST_prot_pool_exchange"]

print(f"Loading {MAT_MODEL} .....", flush=True)
m = load_mat73_cobra(str(MAT_MODEL), var="modelCom")
for solver in ("gurobi", "glpk"):
    try: m.solver = solver; break
    except: continue
print(f"Loaded: {len(m.reactions)} rxns, solver={m.solver.interface.__name__}", flush=True)


PROT_POOL_BOUNDS = [-1000, -100, -50, -20, -10, -5, -2, -1, -0.5, -0.1, -0.05, -0.02, -0.01, -0.005, -0.001]

rows = []
for lb in PROT_POOL_BOUNDS:
    with m:
        for rid in PROT_RXNS:
            if rid in {r.id for r in m.reactions}:
                m.reactions.get_by_id(rid).lower_bound = lb
        m.objective = BIOMASS_ID
        sol = m.optimize()
        mu  = round(sol.objective_value, 8) if sol.status == "optimal" else float("nan")
        rows.append({"prot_pool_lb": lb, "status": sol.status, "mu_h-1": mu})
        print(f"  lb={lb:>8}  {sol.status:<12}  mu={mu:.6f}", flush=True)

df = pd.DataFrame(rows)
print("\n", df.to_string(index=False), flush=True)
RESULTS.mkdir(parents=True, exist_ok=True)
df.to_excel(OUTPUT_XLSX, index=False)
print(f"\nSaved: {OUTPUT_XLSX}")
