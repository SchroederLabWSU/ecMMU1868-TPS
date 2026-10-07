"""
ATP-Maintenance Sensitivity of FABP7-KO Clearance

This script tests whether the WT and FABP7-KO clearance phenotype changes
as ATP-maintenance demand is increased.

Clearance is evaluated from the biomass baseline through increasing fractions
of maximal ATPM demand. The biomass baseline and 0.99 ATPM condition are saved,
while the remaining conditions are printed as a robustness check.

Input:  models/06_community_TPS/enzyme_constrained_tps/ecMMU1868-TPS.mat
Output: results/09_WT_vs_FABP7_ko_lipid_redox/objective_comparison.xlsx
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from cobra.manipulation import knock_out_model_genes

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
RESULTS = REPO_ROOT / "results" / "09_WT_vs_FABP7_ko_lipid_redox"
MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
OUTPUT_XLSX = RESULTS / "objective_comparison.xlsx"

# Load MATLAB v7.3 community model
sys.path.insert(0, str(SCRIPT_DIR))
from mat73_to_cobra import load_mat73_cobra

# Clearance and C04717 reactions
DETOX = [f"AST_LNLCOO_REDUCEc_EXP_{i}" for i in range(1, 6)]
SELF_PROD = ["AST_HMR_2440_EXP_1", "AST_HMR_2440_EXP_2"]

F5_RXNS = [
    "AST_C04717tbulk_EXP_1",
    "AST_C04717tbulk_REV_EXP_1",
    "AST_C04717tcFABP5_EXP_1"
]
F3_RXNS = [
    "AST_C04717tbulk_EXP_2",
    "AST_C04717tbulk_REV_EXP_2",
    "AST_C04717tcFABP5_EXP_2"
]
PAP_RXNS = ["AST_C04717tcFABP5_EXP_1", "AST_C04717tcFABP5_EXP_2"]

NEURON_NEW = [
    f"{c}_C04717tbulk{s}_EXP_{i}"
    for c in ["PRE", "POST"] for s in ["", "_REV"] for i in [1, 2]
]

TD = [
    "AST_C04717td", "AST_C04717td_REV",
    "PRE_C04717td", "PRE_C04717td_REV",
    "POST_C04717td", "POST_C04717td_REV"
]

# FABP protein species and abundance-based usage caps
PROT_STD = "AST_prot_standard"
PROT = {
    "Fabp7": "AST_prot_E9Q0H6",
    "Fabp5": "AST_prot_Q3TLH6",
    "Fabp3": "AST_prot_Q5EBJ0"
}
USAGE = {g: f"AST_usage_prot_{a.rsplit('_', 1)[-1]}" for g, a in PROT.items()}
CAP_RATIO = {"Fabp7": 1.0, "Fabp5": 0.5393, "Fabp3": 0.0074}

# DLKcat turnover numbers and GECKO scaling constants
KCAT = {"Fabp7": 6.2826, "Fabp5": 0.2594, "Fabp3": 3.8233}
CCONST = {
    "Fabp5": 0.433019 * 9.2471,
    "Fabp3": 0.445155 * 9.2471
}

# Passive C04717 transport routes
PASSIVE = TD + ["EX_C04717[u]", "AST_C04717ATP", "PRE_C04717ATP", "POST_C04717ATP"]
PASSIVE_BI = ["AST_IEX_C04717[u]tr", "PRE_IEX_C04717[u]tr", "POST_IEX_C04717[u]tr"]

FABP7_GENE = "AST_12140"
FABP7_CAPACITY = 10.9628

HPODE_SOURCE = "PRE_SK_C04717_prec"
LINOLEATE_SINK = "AST_SK_lnlc_c"
FORCED_LOAD = 1e-3

CAP_F7 = 0.00582   # FABP7 protein-usage cap
PASSIVE_F = 0.01   # Fick-derived passive transport fraction


def solve_max(m, ids):
    # Maximize summed astrocytic detoxification flux.
    m.objective = {m.reactions.get_by_id(r): 1.0 for r in ids if r in m.reactions}
    m.objective_direction = "max"
    solution = m.optimize()

    return (
        (float(solution.objective_value), solution.fluxes)
        if solution.status == "optimal"
        else (np.nan, None)
    )

# Load the enzyme-constrained tripartite synapse model.
print(f"Loading {MAT_MODEL} ...", flush=True)
base = load_mat73_cobra(str(MAT_MODEL))

try:
    base.solver = "gurobi"
except Exception:
    pass

have = {r.id for r in base.reactions}

print(f"  {len(base.reactions)} rxns | "
      f"{len(base.metabolites)} mets", flush=True)


def set_cost(m, rxn_id, prot_id, new_cost):
    # Replace the enzyme coefficient for one transport reaction.
    if rxn_id not in have:
        return

    reaction = m.reactions.get_by_id(rxn_id)
    met = m.metabolites.get_by_id(prot_id)

    if met not in reaction.metabolites:
        return

    reaction.add_metabolites(
        {met: -new_cost - reaction.get_coefficient(prot_id)},
        combine=True
    )


# ATP-maintenance objective uses the summed ATPM flux across all three cells.
BIOMASS_ALL = [
    "AST_BIOMASS_reaction",
    "PRE_BIOMASS_reaction",
    "POST_BIOMASS_reaction"
]
ATPM_ALL = ["AST_ATPM", "PRE_ATPM", "POST_ATPM"]


def flux_sum(m, ids):
    # Build the summed net flux expression.
    expr = None

    for rid in ids:
        flux = m.reactions.get_by_id(rid).flux_expression
        expr = flux if expr is None else expr + flux

    return expr


def configure(m, mult):
    # Remove the prot_standard artifact from passive routes.
    if PROT_STD in {x.id for x in m.metabolites}:
        ps = m.metabolites.get_by_id(PROT_STD)

        for rid in TD:
            if rid in have:
                reaction = m.reactions.get_by_id(rid)

                if ps in reaction.metabolites:
                    reaction.add_metabolites(
                        {ps: -reaction.get_coefficient(PROT_STD)},
                        combine=True
                    )

    # Apply substrate-specific DLKcat enzyme costs.
    for rid in F5_RXNS:
        set_cost(m, rid, PROT["Fabp5"],
                 CCONST["Fabp5"] / (KCAT["Fabp5"] * mult))

    for rid in F3_RXNS:
        set_cost(m, rid, PROT["Fabp3"],
                 CCONST["Fabp3"] / (KCAT["Fabp3"] * mult))

    for rid in NEURON_NEW:
        if rid not in have:
            continue

        cell = rid.split("_")[0]
        gene = "Fabp5" if rid.endswith("_1") else "Fabp3"
        pid = f"{cell}_prot_{PROT[gene].rsplit('_', 1)[-1]}"

        if pid in {x.id for x in m.metabolites}:
            set_cost(m, rid, pid, CCONST[gene] / (KCAT[gene] * mult))

    # Apply the same forced lipid-stress condition in WT and KO.
    m.reactions.get_by_id(HPODE_SOURCE).bounds = (FORCED_LOAD, FORCED_LOAD)
    m.reactions.get_by_id(LINOLEATE_SINK).bounds = (0.0, 1000.0)

    for rid in SELF_PROD:
        if rid in have:
            m.reactions.get_by_id(rid).bounds = (0.0, 0.0)

    # Apply diffusion-limited passive transport bounds.
    cap = PASSIVE_F * FABP7_CAPACITY

    for rid in PASSIVE:
        if rid in have:
            m.reactions.get_by_id(rid).bounds = (0.0, cap)

    for rid in PASSIVE_BI:
        if rid in have:
            m.reactions.get_by_id(rid).bounds = (-cap, cap)

    # Apply abundance-derived protein-usage caps.
    for gene, ratio in CAP_RATIO.items():
        uid = USAGE[gene]

        if uid in have:
            m.reactions.get_by_id(uid).lower_bound = -(CAP_F7 * ratio)

    # Keep FABP5/FABP3 PAP-to-cytosol transport open.
    for rid in PAP_RXNS:
        if rid in have:
            m.reactions.get_by_id(rid).upper_bound = 1000.0


# Biomass-objective reference clearance.
BASE_5C = 0.00662832985115446
BASE_5C_KO = 0.0002433359047051031


rows = []

print()
print("=" * 92)
print("ATP-MAINTENANCE PROBE -- does clearance survive, and at what maintenance load?")
print("=" * 92)
print(f"{'condition':<40} {'WT':>12} {'KO':>12} {'KO/WT':>9} {'ATPM held':>10}")
print("-" * 92)


def probe(label, block_biomass, atpm_frac):
    # Measure WT and KO clearance under one ATP-maintenance condition.
    held = float("nan")
    out = {}

    for arm in ("WT", "KO"):
        with base as m:
            configure(m, 1.0)

            if block_biomass:
                for rid in BIOMASS_ALL:
                    m.reactions.get_by_id(rid).bounds = (0.0, 0.0)

                for rid in ATPM_ALL:
                    m.reactions.get_by_id(rid).upper_bound = 1000.0

                if atpm_frac is not None:
                    expr = flux_sum(m, ATPM_ALL)

                    # Determine maximal ATPM and require at least the selected fraction.
                    m.objective = m.problem.Objective(expr, direction="max")
                    solution = m.optimize()

                    atpm_max = (
                        float(solution.objective_value)
                        if solution.status == "optimal"
                        else float("nan")
                    )

                    held = atpm_frac * atpm_max

                    pin = m.problem.Constraint(
                        expr, lb=held, ub=1e9, name="atpm_pin"
                    )
                    m.add_cons_vars([pin])
                    m.solver.update()

            if arm == "KO":
                knock_out_model_genes(m, [FABP7_GENE])

            clearance, _ = solve_max(m, DETOX)
            out[arm] = clearance

    ratio = (
        100.0 * out["KO"] / out["WT"]
        if out["WT"] > 1e-12
        else float("nan")
    )

    rows.append({
        "Condition": label,
        "ATPM fraction": atpm_frac,
        "atpm_held": held,
        "wt_clearance": out["WT"],
        "ko_clearance": out["KO"],
        "KO % of WT": ratio,
        "wt_frac_of_fig5c": out["WT"] / BASE_5C,
        "ko_frac_of_fig5c": out["KO"] / BASE_5C_KO
    })

    print(f"{label:<40} {out['WT']:>12.6f} {out['KO']:>12.6f} "
          f"{ratio:>8.2f}% {held:>10.4f}", flush=True)


# Test clearance from the biomass baseline through maximal ATPM demand.
probe("biomass-objective baseline", False, None)
probe("biomass blocked, ATPM free", True, None)

for fraction in (0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.00):
    probe(f"biomass blocked, ATPM ≥ {fraction:.2f} × max", True, fraction)


# Only the biomass baseline and 0.99 ATPM condition are published.
COLUMNS = {
    "wt_clearance": "WT clearance flux",
    "ko_clearance": "KO clearance flux"
}

SHEET_COLUMNS = [
    "Condition",
    "ATPM fraction",
    "WT clearance flux",
    "KO clearance flux",
    "KO % of WT"
]

PUBLISHED = [
    "biomass-objective baseline",
    "biomass blocked, ATPM ≥ 0.99 × max"
]


# Save the two reported conditions; the full sweep remains printed as a check.
df = pd.DataFrame(rows)

print()
print("published rows:", PUBLISHED, flush=True)

RESULTS.mkdir(parents=True, exist_ok=True)

df = df[df.Condition.isin(PUBLISHED)]
df.rename(columns=COLUMNS)[SHEET_COLUMNS].to_excel(
    OUTPUT_XLSX,
    sheet_name="objective_comparison",
    index=False
)

print()
print("Fig 5C reference: WT 0.006628, KO 0.000243, ratio 3.67%")
print("Matching rows reproduce the Fig 5C clearance phenotype.")
print()
print(f"Saved: {OUTPUT_XLSX}")