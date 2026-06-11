"""
 This quantifies FABP7-dependent clearance of neuron-derived peroxidized lipid
(C04717) in the enzyme-constrained tripartite synapse model.

A defined peroxidized-lipid load is delivered to the astrocytic process, and
clearance is measured as the FABP7-mediated import flux (AST_C04717tc).

WT: clearance increases with load until reaching the maximal FABP7 capacity.
KO: clearance is abolished because FABP7-mediated import is blocked.

Output:
results/09_WT_vs_FABP7_ko_lipid_redox4/FABP7_clearance_capacity.xlsx
"""

from pathlib import Path
import pandas as pd
import numpy as np
import scipy.sparse as sp
import h5py
import gurobipy as gp
from gurobipy import GRB

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
RESULTS = REPO_ROOT / "results" / "09_WT_vs_FABP7_ko_lipid_redox4"
OUT_XLSX = RESULTS / "FABP7_clearance_capacity.xlsx"

FABP7_IMPORT_RXN = "AST_C04717tc"  # FABP7 import [pap]->[c]; objective; deleted in KO
PAP_POOL_MET = "AST_C04717[pap]"  # astrocytic-process pool where the load is delivered
CLEFT_TO_PAP_RXN = "AST_C04717tpap"  # cleft->PAP transport; blocked so the load is the only supply
LINOLEATE_SINK = "AST_SK_lnlc_c"   # drain for repaired linoleate
# FABP7-dependent reactions removed in the knockout (gene 12140 + dependent bypass routes)
KO_RXNS = ["AST_C04717tc", "AST_C04717td", "AST_C04717td_REV", "AST_C04717ATP"]
LOAD_GRID = [0.1, 0.3, 0.5, 1.0, 3.0, 5.0, 10.0, 11.0, 12.0, 15.0, 30.0, 100.0]


def load_matrices(path):
    f = h5py.File(str(path), "r")
    model = f["modelCom"]
    def deref(ref):
        return "".join(chr(c) for c in np.array(f[ref]).ravel())
    rxn_names = [deref(r) for r in np.array(model["rxns"]).ravel()]
    met_names = [deref(r) for r in np.array(model["mets"]).ravel()]
    lower_bounds = np.array(model["lb"]).ravel().astype(float)
    upper_bounds = np.array(model["ub"]).ravel().astype(float)
    rhs = np.array(model["b"]).ravel().astype(float)
    S_group = model["S"]
    S = sp.csc_matrix((np.array(S_group["data"]).ravel(),
            np.array(S_group["ir"]).ravel().astype(int),
            np.array(S_group["jc"]).ravel().astype(int)),
            shape=(len(met_names), len(rxn_names)))
    return rxn_names, met_names, lower_bounds, upper_bounds, rhs, S


def main():
    print(f"Loading {MAT_MODEL.name} ...", flush=True)
    rxn_names, met_names, lower_bounds, upper_bounds, rhs, S = load_matrices(MAT_MODEL)
    rxn_index = {name: i for i, name in enumerate(rxn_names)}
    met_index = {name: i for i, name in enumerate(met_names)}
    print(f"  {len(rxn_names)} reactions, {len(met_names)} metabolites", flush=True)

    import_col = rxn_index[FABP7_IMPORT_RXN]
    pap_row = met_index[PAP_POOL_MET]
    # augment S with a metered SOURCE and a BURDEN sink 
    source_vec = sp.lil_matrix((S.shape[0], 1)); source_vec[pap_row, 0] = 1.0
    burden_vec = sp.lil_matrix((S.shape[0], 1)); burden_vec[pap_row, 0] = -1.0
    S_aug = sp.hstack([S, source_vec, burden_vec]).tocsr()
    burden_col = len(rxn_names) + 1    # columns appended in order: [source, burden]
    n_cols = S_aug.shape[1]

    def clearance(load, knockout):
        lb = np.append(lower_bounds.copy(), [load, 0.0])     # source forced to load; burden free
        ub = np.append(upper_bounds.copy(), [load, 1000.0])
        ub[rxn_index[LINOLEATE_SINK]] = 1000.0
        lb[rxn_index[CLEFT_TO_PAP_RXN]] = 0.0; ub[rxn_index[CLEFT_TO_PAP_RXN]] = 0.0
        if knockout:
            for r in KO_RXNS:
                lb[rxn_index[r]] = 0.0; ub[rxn_index[r]] = 0.0
        model = gp.Model(); model.Params.OutputFlag = 0
        flux = model.addMVar(n_cols, lb=lb, ub=ub)
        model.addConstr(S_aug @ flux == rhs)
        model.setObjective(flux[import_col], GRB.MAXIMIZE)   # maximize FABP7-mediated clearance
        model.optimize()
        if model.Status != GRB.OPTIMAL:
            return float("nan"), float("nan")
        return flux.X[import_col], flux.X[burden_col]

    rows = []
    print(f"\n{'load':>7} | {'WT_clear':>9} {'WT_burden':>9} | {'KO_clear':>9} {'KO_burden':>9}", flush=True)
    for load in LOAD_GRID:
        wt_clear, wt_burden = clearance(load, knockout=False)
        ko_clear, ko_burden = clearance(load, knockout=True)
        rows.append((load, wt_clear, wt_burden, ko_clear, ko_burden))
        print(f"{load:>7} | {wt_clear:>9.4f} {wt_burden:>9.4f} | {ko_clear:>9.4f} {ko_burden:>9.4f}", flush=True)

    capacity = max(row[1] for row in rows if not np.isnan(row[1]))
    print(f"\nWT FABP7 clearance capacity (plateau) = {capacity:.4f} mmol/gDW/h", flush=True)

    RESULTS.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=["load_C04717", "WT_clearance", "WT_burden", "KO_clearance", "KO_burden"])
    df.to_excel(OUT_XLSX, index=False)
    print(f"Saved: {OUT_XLSX}", flush=True)

if __name__ == "__main__":
    main()
