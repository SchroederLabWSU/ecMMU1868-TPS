"""
Protein-Pool Sensitivity Analysis of Lactate-Coupling Onset in ecMMU1868-TPS

The three GECKO protein-pool bounds are scaled together across 25 values from
0.1 to 22 mmol/gDW/h (0.01x to 2.2x the baseline bound of 10).

At each bound, maximum community growth and the lactate-coupling onset are
recomputed under the biomass objective.

Input:  models/06_community_TPS/enzyme_constrained_tps/ecMMU1868-TPS.mat
Output: results/08_lactate_envelopes_ANLS/03_protein_pool_sweep/protein_pool_sweep.xlsx
"""

import sys
import time
import numpy as np
import pandas as pd
from pathlib import Path

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
RESULTS = REPO_ROOT / "results" / "08_lactate_envelopes_ANLS" / "03_protein_pool_sweep"
MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
OUTPUT_XLSX = RESULTS / "protein_pool_sweep.xlsx"

# Load MATLAB v7.3 community model
sys.path.insert(0, str(SCRIPT_DIR.parent))
from mat73_to_cobra import load_mat73_cobra

BIOMASS = "AST_BIOMASS_reaction"
LACTATE = "EX_lac__L[u]"
POOLS = [
    "AST_prot_pool_exchange",
    "PRE_prot_pool_exchange",
    "POST_prot_pool_exchange"
]

# Protein-pool sweep
BASELINE_POOL = 10.0
GRID = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1, 1.5, 2, 3, 5, 7, 9, 10, 12, 14, 15,
    16, 16.5, 17, 18, 19, 20, 21, 22
]

# Solver and coupling settings
SOLVER_TOLERANCE = 1e-7
REL_TOL = 0.01                # min lactate counts as obligate above 1% of max
BISECT_STEPS = 9
MU_FLOOR = 1e-6               # below this the onset bisection is degenerate
TOP_BACKOFF = [1.0, 0.999, 0.998, 0.995, 0.99]

# Stoichiometric-model growth ceiling
STOICH_CEILING = 0.104826
CEILING_TOL = 1e-4
MONO_RTOL = 1e-6               # Small slack so ceiling points are not falsely flagged


def solve(m, rxn, direction):
    # Optimize one reaction and return NaN if the solve is not optimal.
    m.objective = {m.reactions.get_by_id(rxn): 1.0}
    m.objective_direction = direction
    solution = m.optimize()

    return (float(solution.objective_value)
            if solution.status == "optimal" else float("nan"))


def coupling_onset(m, max_growth):
    # Find the lowest growth fraction where lactate becomes obligatory.
    biomass = m.reactions.get_by_id(BIOMASS)

    def probe(frac):
        # Test whether minimum lactate exceeds 1% of the feasible maximum.
        growth = frac * max_growth
        biomass.bounds = (growth - 1e-9, growth + 1e-9)

        high = solve(m, LACTATE, "max")
        low = solve(m, LACTATE, "min")

        if np.isnan(high) or np.isnan(low):
            return None
        if high <= 1e-9:
            return False

        return low > REL_TOL * high

    # Exact maximum growth can be numerically fragile, so back off if needed.
    top_frac, top_ok = None, False

    for frac in TOP_BACKOFF:
        ok = probe(frac)

        if ok is None:
            continue

        top_frac, top_ok = frac, ok
        break

    if top_frac is None:
        return float("nan"), "no probe near the top would solve"

    if not top_ok:
        return float("nan"), "no coupling: minimum lactate stays near zero"

    # Locate the coupling onset by bisection.
    low_frac, high_frac = 0.0, top_frac

    for _ in range(BISECT_STEPS):
        mid = 0.5 * (low_frac + high_frac)
        ok = probe(mid)

        if ok:
            high_frac = mid
        else:
            low_frac = mid

    return high_frac, "ok"


# Load the enzyme-constrained tripartite synapse model.
t0 = time.time()
print(f"Loading {MAT_MODEL} ...", flush=True)

model = load_mat73_cobra(str(MAT_MODEL), var="modelCom")

try:
    model.solver = "gurobi"
except Exception:
    pass

model.tolerance = SOLVER_TOLERANCE

# Confirm that all three cell-specific protein pools are present.
missing = [p for p in POOLS if p not in {r.id for r in model.reactions}]
assert not missing, f"pool reactions absent from the model: {missing}"

print(f"  {len(model.reactions)} rxns, "
      f"solver={model.solver.interface.__name__}", flush=True)


# Sweep all three protein-pool bounds together.
rows = []

for bound in GRID:
    t_point = time.time()

    # Each bound starts from the same unmodified model.
    with model as m:
        for rxn_id in POOLS:
            m.reactions.get_by_id(rxn_id).bounds = (-float(bound), 0.0)

        m.reactions.get_by_id(BIOMASS).bounds = (0.0, 1000.0)

        max_growth = solve(m, BIOMASS, "max")
        at_ceiling = abs(max_growth - STOICH_CEILING) < CEILING_TOL

        if np.isnan(max_growth) or max_growth < MU_FLOOR:
            onset = float("nan")
            note = (
                "growth infeasible at this pool bound"
                if np.isnan(max_growth)
                else f"max growth below {MU_FLOOR:g}; onset undefined"
            )
        else:
            onset, note = coupling_onset(m, max_growth)

    rows.append({
        "Protein pool bound": float(bound),
        "Pool relative to baseline": float(bound) / BASELINE_POOL,
        "Max community growth": max_growth,
        "Enzyme constraint binding": not at_ceiling,
        "Coupling onset fraction": onset,
        "Obligate lactate coupling": bool(not np.isnan(onset)),
        "Growth at onset": onset * max_growth if not np.isnan(onset) else float("nan"),
    })

    shown = "none" if np.isnan(onset) else f"{onset:.4f}"

    print(f"  bound={bound:>5g}  max growth={max_growth:.6f}"
          f"{'  [stoichiometric ceiling]' if at_ceiling else ''}"
          f"  onset={shown}  [{note}]  {time.time() - t_point:.1f}s",
          flush=True)


sweep = (
    pd.DataFrame(rows)
    .sort_values("Protein pool bound")
    .reset_index(drop=True)
)


# Maximum growth should not decrease as the protein budget increases.
tail_min = sweep["Max community growth"][::-1].cummin()[::-1].shift(-1)
bad = (sweep["Max community growth"] > tail_min * (1.0 + MONO_RTOL)).fillna(False)

assert not bad.any(), (
    "growth monotonicity check failed; the LP did not converge at "
    + ", ".join(f"{v:g}" for v in sweep.loc[bad, "Protein pool bound"])
)

print("\nall bounds pass the growth monotonicity check", flush=True)


# Use "NA" to distinguish absence of coupling in the exported table.
for column in ("Coupling onset fraction", "Growth at onset"):
    sweep[column] = sweep[column].astype(object).where(
        sweep[column].notna(), "NA"
    )


# Save sweep results.
RESULTS.mkdir(parents=True, exist_ok=True)
sweep.to_excel(OUTPUT_XLSX, sheet_name="01_sweep", index=False)

print(f"total runtime {time.time() - t0:.0f}s", flush=True)
print(f"\nSaved: {OUTPUT_XLSX}\n")
print(sweep.to_string(index=False))