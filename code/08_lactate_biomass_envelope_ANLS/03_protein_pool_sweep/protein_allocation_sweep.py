"""
Enzyme Allocation Across Protein-Pool Regimes in ecMMU1868-TPS

This script tracks enzyme allocation across the protein-pool regimes and
measures how protein-budget composition changes between bounds.

Per-enzyme usage is taken from the GECKO usage reactions: <CELL>_usage_prot_<UniProt>.

Input:  models/06_community_TPS/enzyme_constrained_tps/ecMMU1868-TPS.mat
Output: results/08_lactate_envelopes_ANLS/03_protein_pool_sweep/protein_allocation_sweep.xlsx
"""

import sys
import time
import numpy as np
import pandas as pd
import h5py
from pathlib import Path
from cobra.flux_analysis import pfba

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
RESULTS = REPO_ROOT / "results" / "08_lactate_envelopes_ANLS" / "03_protein_pool_sweep"
MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
OUTPUT_XLSX = RESULTS / "protein_allocation_sweep.xlsx"

# Load MATLAB v7.3 community model
sys.path.insert(0, str(SCRIPT_DIR.parent))
from mat73_to_cobra import load_mat73_cobra, _read_cell_of_str, _decode_str

CELLS = ["AST", "PRE", "POST"]
BIOMASS = "AST_BIOMASS_reaction"
POOLS = [f"{c}_prot_pool_exchange" for c in CELLS]

# Match only GECKO protein-usage reactions.
USAGE_TOKEN = "_usage_prot_"

# Bounds included in regimes II-IV of the protein-pool sweep.
POOL_BOUNDS = [
    1.5, 2, 3, 5,               # regime II
    7, 9, 10, 12, 14, 15, 16,   # regime III
    16.5, 17, 18                # regime IV
]

REGIME_EDGES = [1.25, 7.0, 16.25, 18.5]
REPRESENTATIVE = {"II": 3.0, "III": 12.0, "IV": 16.5}
COMPARISONS = [("II", "III"), ("III", "IV"), ("II", "IV")]

# Solver and reporting settings
SOLVER_TOLERANCE = 1e-7
USAGE_FLOOR = 1e-9
MOVER_FLOOR = 1e-6
PFBA_FRACTIONS = [1.0, 0.9999, 0.999] # relax the growth pin if pFBA fails


def regime_of(bound):
    # Assign each protein-pool bound to its regime.
    for edge, name in zip(REGIME_EDGES, ["I", "II", "III", "IV"]):
        if bound < edge:
            return name
    return "V"


# Load the enzyme-constrained tripartite synapse model.
t0 = time.time()
print(f"Loading {MAT_MODEL} ...", flush=True)

model = load_mat73_cobra(str(MAT_MODEL), var="modelCom")

try:
    model.solver = "gurobi"
except Exception:
    pass

model.tolerance = SOLVER_TOLERANCE

print(f"  {len(model.reactions)} rxns | {time.time() - t0:.0f}s", flush=True)


# Identify all cell-specific GECKO protein-usage reactions.
usage_rxns = []

for rxn in model.reactions:
    if USAGE_TOKEN not in rxn.id:
        continue

    cell, accession = rxn.id.split(USAGE_TOKEN, 1)
    usage_rxns.append((rxn.id, cell, accession))

print(f"  {len(usage_rxns)} usage reactions "
      f"({len({a for _, _, a in usage_rxns})} accessions x {len(CELLS)} cells)",
      flush=True)

assert len(usage_rxns) > 5000, "usage reactions not found as expected"


# Decode reaction pathways directly from the MATLAB subSystems field.
print("  decoding subSystems ...", flush=True)

pathway_by_rxn = {}

with h5py.File(str(MAT_MODEL), "r") as f:
    group = f["modelCom"]
    rxn_ids = _read_cell_of_str(f, group["rxns"])

    for rxn_id, outer in zip(
        rxn_ids,
        np.asarray(group["subSystems"][()]).flatten()
    ):
        if not outer:
            continue

        try:
            names = [
                _decode_str(f[ref][()])
                for ref in np.asarray(f[outer][()]).flatten()
                if ref
            ]
            names = [n for n in names if n]

            if names:
                pathway_by_rxn[rxn_id] = names[0]

        except Exception:
            continue

print(f"    {len(pathway_by_rxn)} reactions carry a pathway, "
      f"{len(set(pathway_by_rxn.values()))} distinct", flush=True)


# Assign each enzyme to its most common associated pathway.
print("  mapping enzymes to pathways ...", flush=True)

enzyme_meta = {}

for rxn_id, cell, accession in usage_rxns:
    pathway = ""

    try:
        met = model.metabolites.get_by_id(f"{cell}_prot_{accession}")
        names = [
            pathway_by_rxn.get(r.id, "")
            for r in met.reactions
            if USAGE_TOKEN not in r.id
        ]
        names = [n for n in names if n]

        if names:
            pathway = pd.Series(names).value_counts().index[0]

    except KeyError:
        pass

    enzyme_meta[rxn_id] = {
        "cell": cell,
        "accession": accession,
        "Pathway": pathway
    }

print(f"    {sum(1 for v in enzyme_meta.values() if v['Pathway'])} of "
      f"{len(enzyme_meta)} enzymes mapped", flush=True)


# Calculate enzyme allocation across protein-pool bounds.
allocation_rows, summary_rows = [], []

for bound in POOL_BOUNDS:
    regime = regime_of(bound)

    with model as m:
        for rxn_id in POOLS:
            m.reactions.get_by_id(rxn_id).bounds = (-float(bound), 0.0)

        m.reactions.get_by_id(BIOMASS).bounds = (0.0, 1000.0)
        m.objective = {m.reactions.get_by_id(BIOMASS): 1.0}
        m.objective_direction = "max"

        solution = m.optimize()
        assert solution.status == "optimal", f"bound {bound}: {solution.status}"

        max_growth = float(solution.objective_value)

        # Relax the growth fraction slightly if pFBA fails at the exact optimum.
        fluxes = None

        for fraction in PFBA_FRACTIONS:
            try:
                fluxes = pfba(m, fraction_of_optimum=fraction).fluxes
                break
            except Exception:
                continue

        assert fluxes is not None, f"bound {bound}: no pFBA solution"

        active = {}

        for cell in CELLS:
            draw = -float(fluxes.get(f"{cell}_prot_pool_exchange", 0.0))
            total, n_active = 0.0, 0

            for rxn_id, meta in enzyme_meta.items():
                if meta["cell"] != cell:
                    continue

                usage = -float(fluxes.get(rxn_id, 0.0))
                total += usage

                if usage <= USAGE_FLOOR:
                    continue

                n_active += 1

                allocation_rows.append({
                    "Protein pool bound": bound,
                    "Regime": regime,
                    "Cell type": cell,
                    "Protein accession": meta["accession"],
                    "Pathway": meta["Pathway"],
                    "Protein usage": usage,
                    "Protein pool share": usage / draw,
                })

            active[cell] = n_active

            # Protein usages should sum to the cell's pool draw.
            assert abs(total - draw) < 1e-6, \
                f"bound {bound}, {cell}: usage sum {total} != draw {draw}"

            # Pool draw cannot exceed the imposed bound.
            assert draw <= bound * (1.0 + 1e-6), \
                f"bound {bound}, {cell}: draw {draw} exceeds the bound"

    summary_rows.append({
        "Protein pool bound": bound,
        "Regime": regime,
        "Enzymes active AST": active["AST"],
        "Enzymes active PRE": active["PRE"],
        "Enzymes active POST": active["POST"],
    })

    print(f"  bound={bound:>5g}  regime {regime:<3}  "
          f"max growth={max_growth:.6f}  active "
          + " ".join(f"{c} {active[c]}" for c in CELLS), flush=True)


summary = pd.DataFrame(summary_rows)
allocation = pd.DataFrame(allocation_rows)


# Measure protein-budget reallocation between consecutive bounds.
change_rows = []

for cell in CELLS:
    wide = (
        allocation[allocation["Cell type"] == cell]
        .pivot_table(
            index="Protein accession",
            columns="Protein pool bound",
            values="Protein pool share",
            aggfunc="sum"
        )
        .fillna(0.0)
    )

    bounds = sorted(wide.columns)

    for previous, current in zip(bounds, bounds[1:]):
        before, after = regime_of(previous), regime_of(current)

        change_rows.append({
            "From protein pool bound": previous,
            "To protein pool bound": current,
            "Cell type": cell,
            "Percent budget reallocated":
                100.0 * 0.5 * float(np.abs(wide[current] - wide[previous]).sum()),
            "Regime transition":
                before if before == after else f"{before} - {after}",
        })

change = pd.DataFrame(change_rows)


# Aggregate enzyme allocation by pathway.
pathway = (
    allocation[allocation["Pathway"] != ""]
    .groupby(
        ["Protein pool bound", "Regime", "Cell type", "Pathway"],
        as_index=False
    )
    .agg(**{
        "Protein usage": ("Protein usage", "sum"),
        "Protein pool share": ("Protein pool share", "sum")
    })
    .sort_values(
        ["Cell type", "Protein pool bound", "Protein pool share"],
        ascending=[True, True, False]
    )
)


# Compare enzyme shares between representative regime bounds.
mover_rows = []

for cell in CELLS:
    matrix = (
        allocation[allocation["Cell type"] == cell]
        .pivot_table(
            index=["Protein accession", "Pathway"],
            columns="Protein pool bound",
            values="Protein pool share",
            aggfunc="sum"
        )
        .fillna(0.0)
    )

    for first, second in COMPARISONS:
        a, b = REPRESENTATIVE[first], REPRESENTATIVE[second]

        for key, delta in (matrix[b] - matrix[a]).items():
            if abs(delta) < MOVER_FLOOR:
                continue

            mover_rows.append({
                "Cell type": cell,
                "Regime transition": f"{first} - {second}",
                "Protein accession": key[0],
                "Pathway": key[1],
                "Share before": float(matrix.loc[key, a]),
                "Share after": float(matrix.loc[key, b]),
                "Change in share": float(delta),
            })

movers = pd.DataFrame(mover_rows).sort_values(
    ["Cell type", "Regime transition", "Change in share"],
    key=lambda s: s.abs() if s.name == "Change in share" else s,
    ascending=[True, True, False]
)


# Save allocation and regime-comparison results.
RESULTS.mkdir(parents=True, exist_ok=True)

with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as wb:
    summary.to_excel(wb, sheet_name="01_summary", index=False)
    change.to_excel(wb, sheet_name="02_composition_change", index=False)
    allocation.to_excel(wb, sheet_name="03_allocation", index=False)
    pathway.to_excel(wb, sheet_name="04_pathway", index=False)
    movers.to_excel(wb, sheet_name="05_enzyme_share_changes", index=False)

print(f"\nSaved: {OUTPUT_XLSX}  ({time.time() - t0:.0f}s)\n")
print(summary.to_string(index=False))