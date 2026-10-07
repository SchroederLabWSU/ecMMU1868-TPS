"""
Lactate Production Envelopes for iMMU1868-TPS (stoichiometric community)

This script computes the minimum and maximum feasible lactate exchange rates
across the full range of two objectives: community growth and ATP maintenance.

You need to copy the .mat model from Zenodo to run this locally,
since the model file is too large for GitHub. Or run the .m files in
code/07_steadyCom_community/01_stoichiometric_TPS_community to generate it locally.

Input:  models/06_community_TPS/stoichiometric_tps/iMMU1868-TPS.mat
Output: results/08_lactate_envelopes_ANLS/01_envelopes/lactate_envelopes_iMMU1868-TPS.xlsx
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]
RESULTS = REPO_ROOT / "results" / "08_lactate_envelopes_ANLS" / "01_envelopes"
MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "stoichiometric_tps" / "iMMU1868-TPS.mat"
OUTPUT_XLSX = RESULTS / "lactate_envelopes_iMMU1868-TPS.xlsx"

# Load MATLAB v7.3 community model
sys.path.insert(0, str(SCRIPT_DIR.parent))
from mat73_to_cobra import load_mat73_cobra

# Reactions used for the two objectives and lactate envelope
BIOMASS = "AST_BIOMASS_reaction"
BIOMASS_ALL = ["AST_BIOMASS_reaction", "PRE_BIOMASS_reaction", "POST_BIOMASS_reaction"]
ATPM_ALL = ["AST_ATPM", "PRE_ATPM", "POST_ATPM"]
LACTATE = "EX_lac__L[u]"

# Numerical tolerance used when fixing the objective level
TOL = 1e-6
PIN_CEILING = 1e6

# Output column names for each objective
COLUMNS = {
    "biomass": ("Fraction of max biomass", "Biomass flux"),
    "ATPM": ("Fraction of max ATPM", "ATPM flux"),
}


def fraction_grid():
    # Finer spacing is used near maximal demand where coupling may begin.
    return np.unique(np.concatenate([
        np.arange(0.0, 0.60, 0.05),
        np.arange(0.60, 0.80, 0.02),
        np.arange(0.80, 1.0001, 0.005),
    ]))


def flux_sum(model, reaction_ids):
    # Build the summed net flux expression for several reactions.
    expr = 0
    for rid in reaction_ids:
        r = model.reactions.get_by_id(rid)
        expr = expr + r.forward_variable - r.reverse_variable
    return expr


def solve(model, objective, direction):
    # Optimize one objective and return NaN if the solve is not optimal.
    model.objective = objective
    model.objective_direction = direction
    solution = model.optimize()

    return (float(solution.objective_value)
            if solution.status == "optimal" else float("nan"))


def envelope(model, objective):
    # Compute the minimum and maximum feasible lactate flux
    # across the selected biomass or ATPM objective range.
    lactate = model.reactions.get_by_id(LACTATE)
    fraction_column, flux_column = COLUMNS[objective]
    pin = None

    if objective == "biomass":
        # Biomass objective used for the growth-based envelope.
        model.reactions.get_by_id(BIOMASS).bounds = (0.0, 1000.0)
        target = model.reactions.get_by_id(BIOMASS)

        obj_min = 0.0
        obj_max = solve(model, target, "max")
        label, unit = "growth", "h⁻¹"

    else:
        # ATPM objective: block biomass production in all three cell types.
        for rid in BIOMASS_ALL:
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)

        # Release ATPM upper bounds before determining the feasible range.
        for rid in ATPM_ALL:
            model.reactions.get_by_id(rid).upper_bound = 1000.0

        target = flux_sum(model, ATPM_ALL)
        obj_max = solve(model, target, "max")
        obj_min = solve(model, target, "min")

        # ATPM is the sum of three reactions, so it is fixed with a constraint.
        pin = model.problem.Constraint(
            target,
            lb=obj_min,
            ub=PIN_CEILING,
            name="atpm_pin"
        )
        model.add_cons_vars([pin])
        model.solver.update()

        label, unit = "ATP maintenance", "mmol gDW⁻¹ h⁻¹"

    print(f"  {label} range: {obj_min:.6f} to {obj_max:.6f} {unit}", flush=True)

    rows = []

    # Sweep from minimum to maximum objective demand.
    for i, fraction in enumerate(fraction_grid(), 1):
        level = obj_min + fraction * (obj_max - obj_min)

        if pin is None:
            # Biomass is a single reaction, so its bounds can be set directly.
            target.bounds = (level - TOL, level + TOL)

        else:
            # Reset the upper bound first so lb never exceeds ub during the update.
            pin.ub = PIN_CEILING
            pin.lb = level - TOL
            pin.ub = level + TOL

        # Lactate envelope at the current objective level.
        lac_max = solve(model, lactate, "max")
        lac_min = solve(model, lactate, "min")

        rows.append({
            fraction_column: round(float(fraction), 6),
            flux_column: level,
            "Min lactate flux": lac_min,
            "Max lactate flux": lac_max
        })

        if i % 10 == 0 or i == 1:
            print(f"  [{i:3d}] frac={fraction:.3f}  obj={level:.4f}  "
                  f"min_lac={lac_min:.4f}  max_lac={lac_max:.4f}", flush=True)

    # Remove the temporary ATPM constraint after the sweep.
    if pin is not None:
        model.remove_cons_vars([pin])
        model.solver.update()

    return pd.DataFrame(rows)


# Load the stoichiometric tripartite synapse model.
print(f"Loading {MAT_MODEL} ...", flush=True)
model = load_mat73_cobra(str(MAT_MODEL), var="modelCom")

# Prefer Gurobi; fall back to GLPK if unavailable.
for solver in ("gurobi", "glpk"):
    try:
        model.solver = solver
        break
    except:
        continue

print(f"  {len(model.reactions)} rxns, {len(model.metabolites)} mets, "
      f"solver={model.solver.interface.__name__}", flush=True)

# Run each objective from the same unmodified starting model.
envelopes = {}

for objective in COLUMNS:
    print(f"\n{objective} objective", flush=True)

    with model as m:
        envelopes[objective] = envelope(m, objective)

# Save biomass and ATPM envelopes to separate sheets.
RESULTS.mkdir(parents=True, exist_ok=True)

with pd.ExcelWriter(OUTPUT_XLSX) as writer:
    for objective, df in envelopes.items():
        df.to_excel(writer, sheet_name=objective, index=False)

print(f"\nSaved: {OUTPUT_XLSX}")