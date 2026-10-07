"""
kcat Bootstrap Analysis of Lactate-Coupling Onset in ecMMU1868-TPS

Reactions with both BRENDA and DLKcat turnover numbers are randomly assigned
one of the two values across 100 parameterizations. Lactate-coupling onset is
then calculated under either the biomass or ATP-maintenance objective.

The same random seeds are used for both objectives so results can be paired.

Run with:
    --objective biomass
    --objective atpm

Input:  models/06_community_TPS/enzyme_constrained_tps/ecMMU1868-TPS.mat
        results/08_lactate_envelopes_ANLS/02_kcat_bootstrap/kcat_parameter_sheet.xlsx
Output: results/08_lactate_envelopes_ANLS/02_kcat_bootstrap/kcat_bootstrap_<objective>.xlsx
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
RESULTS = REPO_ROOT / "results" / "08_lactate_envelopes_ANLS" / "02_kcat_bootstrap"
MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
PARAMS = RESULTS / "kcat_parameter_sheet.xlsx"

# Allow a local model path when needed.
MAT_MODEL = Path(os.environ.get("ECMMU_MAT", MAT_MODEL))

def output_path(objective):
    return RESULTS / f"kcat_bootstrap_{objective}.xlsx"

# Load MATLAB v7.3 community model
sys.path.insert(0, str(SCRIPT_DIR.parent))
from mat73_to_cobra import load_mat73_cobra

BIOMASS_ID = "AST_BIOMASS_reaction"
LACTATE_ID = "EX_lac__L[u]"
ATPM_RXNS = ["AST_ATPM", "PRE_ATPM", "POST_ATPM"]

# Bootstrap and coupling settings
N_RUNS = 100
BASE_SEED = 20260818
TOL = 1e-7
LACTATE_FRAC = 0.01
BISECTION_STEPS = 9
PIN_NAME = "objective_pin"

# Command-line options
parser = argparse.ArgumentParser()
parser.add_argument("--objective", choices=("biomass", "atpm"), default="biomass")
parser.add_argument("--start", type=int, default=0)
parser.add_argument("--n", type=int, default=N_RUNS)
args = parser.parse_args()


# Load the enzyme-constrained tripartite synapse model.
print(f"Loading {MAT_MODEL} ...", flush=True)
model = load_mat73_cobra(str(MAT_MODEL), var="modelCom")

# Prefer Gurobi; fall back to GLPK if unavailable.
for solver in ("gurobi", "glpk"):
    try:
        model.solver = solver
        break
    except:
        continue

model.solver.configuration.tolerances.feasibility = TOL

print(f"Loaded: {len(model.reactions)} rxns, "
      f"solver={model.solver.interface.__name__}", flush=True)


# Select reactions with both BRENDA and DLKcat turnover numbers.
params = pd.read_excel(PARAMS, sheet_name="01_parameters")

paired = params[
    (params["Source assigned"].astype(str) == "brenda")
    & params["kcat DLKcat cached"].notna()
    & (params["kcat DLKcat cached"] > 0)
    & (params["kcat assigned"] > 0)
].copy()

# Remove pairs with effectively identical turnover numbers.
paired = paired[
    (paired["kcat assigned"] - paired["kcat DLKcat cached"]).abs()
    / paired["kcat assigned"] > 1e-6
]

print(f"Paired kcat reactions available for resampling: {len(paired)}", flush=True)


def objective_expression(m, objective):
    # Biomass uses one reaction; ATPM is summed across the three cell types.
    if objective == "biomass":
        return m.reactions.get_by_id(BIOMASS_ID).flux_expression

    expr = m.reactions.get_by_id(ATPM_RXNS[0]).flux_expression
    for rid in ATPM_RXNS[1:]:
        expr = expr + m.reactions.get_by_id(rid).flux_expression

    return expr


def pin_objective(m, objective, level):
    # Fix the selected objective while optimizing lactate flux.
    existing = m.constraints.get(PIN_NAME) if hasattr(m, "constraints") else None

    if existing is not None:
        m.remove_cons_vars([existing])

    pin = m.problem.Constraint(
        objective_expression(m, objective),
        lb=level, ub=level, name=PIN_NAME
    )
    m.add_cons_vars([pin])
    m.solver.update()


def lactate_range(m, objective, level):
    # Minimum and maximum feasible lactate flux at a fixed objective level.
    pin_objective(m, objective, level)

    m.objective = LACTATE_ID
    m.objective.direction = "min"
    low = m.slim_optimize()

    m.objective.direction = "max"
    high = m.slim_optimize()

    return low, high


def coupling_onset(m, objective, max_level):
    # Find the first objective fraction where lactate becomes obligatory.
    low_frac, high_frac = 0.0, 1.0
    min_lac, max_lac = lactate_range(m, objective, max_level)

    if not (max_lac and min_lac > LACTATE_FRAC * max_lac):
        return None, min_lac, max_lac

    for _ in range(BISECTION_STEPS):
        mid = (low_frac + high_frac) / 2.0
        lo, hi = lactate_range(m, objective, mid * max_level)

        if hi and lo > LACTATE_FRAC * hi:
            high_frac = mid
        else:
            low_frac = mid

    return high_frac, min_lac, max_lac


# Store BRENDA and DLKcat values for each paired reaction.
PAIRS = {
    str(rec["Reaction"]): (
        float(rec["kcat assigned"]),
        float(rec["kcat DLKcat cached"])
    )
    for _, rec in paired.iterrows()
}

# Cache enzyme coefficients for paired reactions present in the model.
COEF, absent, no_prot = {}, [], []

for rxn_id in PAIRS:
    if rxn_id not in model.reactions:
        absent.append(rxn_id)
        continue

    reaction = model.reactions.get_by_id(rxn_id)

    for met in reaction.metabolites:
        if "_prot_" in met.id:
            COEF[rxn_id] = (met.id, reaction.get_coefficient(met.id))
            break
    else:
        no_prot.append(rxn_id)

print(f"Paired in sheet    : {len(PAIRS)}", flush=True)
print(f"Absent from model  : {len(absent)}", flush=True)
print(f"Carrying no _prot_ : {len(no_prot)}", flush=True)
print(f"Resampled          : {len(COEF)}", flush=True)


# Run the turnover-number parameterizations.
rows = []

for run in range(args.start, args.start + args.n):
    seed = BASE_SEED + run
    rng = np.random.default_rng(seed)

    # Each paired reaction has equal probability of using BRENDA or DLKcat.
    take_dlkcat = rng.random(len(COEF)) < 0.5

    with model:
        # Rescale enzyme coefficients only for reactions assigned DLKcat.
        for (rxn_id, (prot_id, original)), use_dl in zip(COEF.items(), take_dlkcat):
            if not use_dl:
                continue

            kcat_assigned, kcat_dlkcat = PAIRS[rxn_id]
            target = original * (kcat_assigned / kcat_dlkcat)

            reaction = model.reactions.get_by_id(rxn_id)
            current = reaction.get_coefficient(prot_id)

            if abs(current - target) > 1e-15:
                reaction.add_metabolites(
                    {model.metabolites.get_by_id(prot_id): target - current},
                    combine=True
                )

        # Set the selected community objective.
        if args.objective == "biomass":
            model.objective = BIOMASS_ID

        else:
            model.reactions.get_by_id(BIOMASS_ID).bounds = (0.0, 0.0)
            model.objective = {model.reactions.get_by_id(r): 1.0 for r in ATPM_RXNS}

        model.objective.direction = "max"
        max_level = model.slim_optimize()

        # Record the ATPM contribution from each cell type.
        cell_atpm = {r: 0.0 for r in ATPM_RXNS}

        if args.objective == "atpm":
            solution = model.optimize()

            if solution.status == "optimal":
                cell_atpm = {r: float(solution.fluxes[r]) for r in ATPM_RXNS}

        onset, min_lac, max_lac = coupling_onset(model, args.objective, max_level)

    # Store results using objective-specific column names.
    if args.objective == "biomass":
        rows.append({
            "Bootstrap run": run,
            "seed": seed,
            "Number of DLKcat reactions": int(take_dlkcat.sum()),
            "Maximum biomass flux": round(max_level, 8),
            "Lactate coupling onset": onset,
            "Lactate min at max biomass": round(min_lac, 6),
            "Lactate max at max biomass": round(max_lac, 6),
            "Couples": onset is not None,
        })

    else:
        rows.append({
            "Bootstrap run": run,
            "seed": seed,
            "Number of DLKcat reactions": int(take_dlkcat.sum()),
            "Max ATPM flux": round(max_level, 8),
            "AST ATPM flux": round(cell_atpm["AST_ATPM"], 6),
            "PRE ATPM flux": round(cell_atpm["PRE_ATPM"], 6),
            "POST ATPM flux": round(cell_atpm["POST_ATPM"], 6),
            "Lactate coupling onset": onset,
            "Lactate min at max ATPM": round(min_lac, 6),
            "Lactate max at max ATPM": round(max_lac, 6),
            "Couples": onset is not None,
        })

    label = f"{onset:.4f}" if onset is not None else "none"
    print(f"  run {run:>3}  max={max_level:.5f}  onset={label}", flush=True)


runs = pd.DataFrame(rows)

# Summarize coupling onset and maximum objective flux.
if args.objective == "atpm":
    flux_column, flux_label = "Max ATPM flux", "median max_ATPM flux"
else:
    flux_column, flux_label = "Maximum biomass flux", "median maximum biomass flux"

onsets = runs["Lactate coupling onset"].dropna()

summary = pd.DataFrame([
    {"statistic": "runs with obligate coupling", "value": int(runs["Couples"].sum())},
    {"statistic": flux_label, "value": round(runs[flux_column].median(), 6)},
    {"statistic": "median lactate-coupling onset", "value": round(onsets.median(), 6)},
    {"statistic": "mean lactate-coupling onset", "value": round(onsets.mean(), 6)},
    {"statistic": "lactate onset SD", "value": round(onsets.std(ddof=1), 6)},
    {"statistic": "5th percentile onset", "value": round(onsets.quantile(0.05), 6)},
    {"statistic": "95th percentile onset", "value": round(onsets.quantile(0.95), 6)},
    {"statistic": "min lactate onset", "value": round(onsets.min(), 6)},
    {"statistic": "max lactate onset", "value": round(onsets.max(), 6)},
])

print("\n", summary.to_string(index=False), flush=True)


# Save the individual runs and summary.
RESULTS.mkdir(parents=True, exist_ok=True)
OUTPUT_XLSX = output_path(args.objective)

with pd.ExcelWriter(OUTPUT_XLSX) as writer:
    runs.to_excel(writer, sheet_name="01_bootstrap_runs", index=False)
    summary.to_excel(writer, sheet_name="02_summary", index=False)


# Pair ATPM and biomass runs using their shared random seeds.
BIOMASS_XLSX = output_path("biomass")

if args.objective == "atpm" and BIOMASS_XLSX.exists():
    biomass_runs = pd.read_excel(BIOMASS_XLSX, sheet_name="01_bootstrap_runs")
    biomass_runs = biomass_runs.rename(columns={"Seed": "seed"})

    paired = runs.merge(
        biomass_runs[["seed", "Lactate coupling onset"]],
        on="seed", how="inner", suffixes=("", "_biomass")
    )

    paired = pd.DataFrame({
        "Bootstrap run": paired["Bootstrap run"],
        "Lactate coupling onset ATPM": paired["Lactate coupling onset"],
        "Lactate coupling onset biomass": paired["Lactate coupling onset biomass"],
    })

    paired["Paired difference"] = (
        paired["Lactate coupling onset ATPM"]
        - paired["Lactate coupling onset biomass"]
    )
    earlier = int((paired["Paired difference"] < 0).sum())

    summary = pd.concat([summary, pd.DataFrame([
        {"statistic": "median biomass-objective onset", "value": round(paired["Lactate coupling onset biomass"].median(), 6)},
        {"statistic": "median paired difference (ATPM-biomass)", "value": round(paired["Paired difference"].median(), 6)},
        {"statistic": "runs where ATPM couples earlier", "value": earlier},
    ])], ignore_index=True)

    # Rewrite the ATPM workbook with the paired comparison included.
    with pd.ExcelWriter(OUTPUT_XLSX) as writer:
        runs.to_excel(writer, sheet_name="01_bootstrap_runs", index=False)
        summary.to_excel(writer, sheet_name="02_summary", index=False)
        paired.to_excel(writer, sheet_name="03_paired", index=False)

    print(f"paired with {len(paired)} biomass runs; "
          f"ATPM couples earlier in {earlier}", flush=True)

print(f"\nSaved: {OUTPUT_XLSX}")