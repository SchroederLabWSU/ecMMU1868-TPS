"""
Model Metadata for Reporter-Metabolite Analysis

Writes the shared model tables used by the reporter-metabolite analysis:
the stoichiometric matrix as metabolite-reaction-coefficient triplets,
reaction metadata, and metabolite metadata.

The .mat is read directly with h5py, so neither MATLAB nor the COBRA Toolbox is
needed to reproduce this step.

Input:  models/06_community_TPS/enzyme_constrained_tps/ecMMU1868-TPS.mat
Output: model_S_triplets.csv, model_rxn_meta.csv, model_met_meta.csv
"""

import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
OUT = REPO_ROOT / "results" / "10_reporter_metabolites"

# Read MATLAB v7.3 fields
sys.path.insert(0, str(SCRIPT_DIR))
from mat73_to_cobra import _read_cell_of_str


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"reading {MODEL.name} ...", flush=True)
    f = h5py.File(MODEL, "r")
    g = f["modelCom"]

    rxns = _read_cell_of_str(f, g["rxns"])
    mets = _read_cell_of_str(f, g["mets"])

    print(f"  {len(rxns)} rxns | {len(mets)} mets", flush=True)

    def cell(name, n):
        # Read a string field and pad with blanks if needed.
        if name not in g:
            print(f"  [warn] {name} absent; writing blanks", flush=True)
            return [""] * n

        values = _read_cell_of_str(f, g[name])
        return values if len(values) == n else (values + [""] * (n - len(values)))[:n]

    # Reaction metadata
    rxn_meta = pd.DataFrame({
        "rxn_id": rxns,
        "rxn_name": cell("rxnNames", len(rxns)),
        "subsystem": cell("subSystems", len(rxns)),
        "grRule": cell("grRules", len(rxns)),
        "lb": np.asarray(g["lb"]).flatten(),
        "ub": np.asarray(g["ub"]).flatten(),
    })

    rxn_meta.to_csv(OUT / "model_rxn_meta.csv", index=False)
    print(f"  wrote model_rxn_meta.csv   ({len(rxn_meta)} rows)", flush=True)

    # Metabolite metadata
    met_meta = pd.DataFrame({
        "met_id": mets,
        "met_name": cell("metNames", len(mets)),
        "formula": cell("metFormulas", len(mets)),
    })

    met_meta.to_csv(OUT / "model_met_meta.csv", index=False)
    print(f"  wrote model_met_meta.csv   ({len(met_meta)} rows)", flush=True)

    # Convert the sparse MATLAB S matrix from CSC format to triplets.
    S = g["S"]
    data = np.asarray(S["data"]).flatten()
    ir = np.asarray(S["ir"]).flatten()
    jc = np.asarray(S["jc"]).flatten()

    rxn_col = np.empty(len(data), dtype=np.int64)

    for j in range(len(jc) - 1):
        rxn_col[jc[j]:jc[j + 1]] = j

    trip = pd.DataFrame({
        "rxn_id": [rxns[j] for j in rxn_col],
        "met_id": [mets[i] for i in ir],
        "coef": data,
    })

    trip.to_csv(OUT / "model_S_triplets.csv", index=False)
    print(f"  wrote model_S_triplets.csv ({len(trip)} nonzeros)", flush=True)

    f.close()
    print("\ndone.")


if __name__ == "__main__":
    main()