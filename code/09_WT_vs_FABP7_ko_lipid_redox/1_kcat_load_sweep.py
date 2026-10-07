"""
FABP Turnover Sensitivity of Peroxidized-Lipid Clearance

Peroxidized lipid (13-HPODE; C04717) is delivered at the presynaptic neuron
and cleared through astrocytic glutathione-dependent reduction.

FABP5 and FABP3 turnover numbers are varied from 0.1x to 10x their DLKcat
predictions. FABP7 is held fixed as the wild-type reference and deleted in
the knockout.

Input:  models/06_community_TPS/enzyme_constrained_tps/ecMMU1868-TPS.mat
Output: results/09_WT_vs_FABP7_ko_lipid_redox/kcat_load_sweep.xlsx
"""

import sys
import math
import numpy as np
import pandas as pd
from pathlib import Path
from cobra import Reaction
from cobra.manipulation import knock_out_model_genes
from openpyxl.styles import Border, Font

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
RESULTS = REPO_ROOT / "results" / "09_WT_vs_FABP7_ko_lipid_redox"
MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
OUTPUT_XLSX = RESULTS / "kcat_load_sweep.xlsx"

# Load MATLAB v7.3 community model
sys.path.insert(0, str(SCRIPT_DIR))
from mat73_to_cobra import load_mat73_cobra

# Clearance and C04717 reactions
DETOX = [f"AST_LNLCOO_REDUCEc_EXP_{i}" for i in range(1, 6)]
FABP7_RXN = "AST_C04717tc"
LINOLEATE_SINK = "AST_SK_lnlc_c"
HPODE_SOURCE = "PRE_SK_C04717_prec"

# FABP5 and FABP3 transport routes
F5 = ["AST_C04717tbulk_EXP_1", "AST_C04717tbulk_REV_EXP_1", "AST_C04717tcFABP5_EXP_1"]
F3 = ["AST_C04717tbulk_EXP_2", "AST_C04717tbulk_REV_EXP_2", "AST_C04717tcFABP5_EXP_2"]
NEURON_NEW = [
    f"{c}_C04717tbulk{s}_EXP_{i}"
    for c in ["PRE", "POST"] for s in ["", "_REV"] for i in [1, 2]
]
NEW_ALL = F5 + F3 + NEURON_NEW

# Passive C04717 transport routes
TD = [
    "AST_C04717td", "AST_C04717td_REV",
    "PRE_C04717td", "PRE_C04717td_REV",
    "POST_C04717td", "POST_C04717td_REV"
]
PASSIVE = TD + ["EX_C04717[u]", "AST_C04717ATP", "PRE_C04717ATP", "POST_C04717ATP"]
PASSIVE_BI = ["AST_IEX_C04717[u]tr", "PRE_IEX_C04717[u]tr", "POST_IEX_C04717[u]tr"]

# FABP protein species and abundance-based usage caps
PROT_STD = "AST_prot_standard"
PROT = {
    "Fabp7": "AST_prot_E9Q0H6",
    "Fabp5": "AST_prot_Q3TLH6",
    "Fabp3": "AST_prot_Q5EBJ0"
}
USAGE = {g: f"AST_usage_prot_{a.rsplit('_', 1)[-1]}" for g, a in PROT.items()}

CAP_RATIO = {"Fabp7": 1.0, "Fabp5": 0.5393, "Fabp3": 0.0074}
CAP_F7 = 0.00582

# DLKcat turnover numbers and original GECKO scaling constants
KCAT_DLK = {"Fabp5": 0.2594, "Fabp3": 3.8233}
CCONST = {
    "Fabp5": 0.433019 * 9.2471,
    "Fabp3": 0.445155 * 9.2471
}

FABP7_GENE = "AST_12140"
FABP7_CAPACITY = 10.9628
PASSIVE_F = 0.01

# Multipliers span ±1 log10 around the DLKcat predictions.
MULTS = [0.1, 0.3162, 1.0, 3.1623, 10.0]
WT_REF_MULT = 1.0

# Delivered-load grid spans the KO and WT clearance plateaus.
LOADS = [
    0.0,
    0.00001, 0.00002, 0.00005, 0.0001,
    0.0002, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005,
    0.006, 0.007, 0.008, 0.010, 0.015, 0.020,
    0.030, 0.045, 0.066, 0.080, 0.100
]


def solve_max(m, ids):
    # Maximize the summed astrocytic detoxification flux.
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
print(f"  {len(base.reactions)} rxns", flush=True)


# Remove the prot_standard artifact from passive diffusion routes.
ps = (
    base.metabolites.get_by_id(PROT_STD)
    if PROT_STD in {m.id for m in base.metabolites}
    else None
)

if ps is not None:
    for rid in TD:
        if rid in have:
            reaction = base.reactions.get_by_id(rid)

            if ps in reaction.metabolites:
                reaction.add_metabolites(
                    {ps: -reaction.get_coefficient(PROT_STD)},
                    combine=True
                )


# Cache the protein species used by the added FABP5/FABP3 routes.
ORIG = {}

for rid in NEW_ALL:
    if rid not in have:
        continue

    reaction = base.reactions.get_by_id(rid)

    for met in reaction.metabolites:
        if "_prot_" in met.id:
            ORIG[rid] = met.id
            break


# FABP7 remains fixed at its baseline model coefficient.
F7_PROT = PROT["Fabp7"]
F7_ORIG = base.reactions.get_by_id(FABP7_RXN).get_coefficient(F7_PROT)
print(f"  FABP7 {FABP7_RXN} enzyme coefficient: {F7_ORIG:.6g}", flush=True)

# Block endogenous HMR_2440 peroxidation so delivered load is controlled directly.
HMR = [r.id for r in base.reactions if "HMR_2440" in r.id]
print(f"  peroxidation reactions to block: {len(HMR)}", flush=True)


def configure(m, mult, load):
    # Rescale FABP5/FABP3 turnover while keeping FABP7 fixed.
    for rid, pid in ORIG.items():
        gene = "Fabp5" if rid.endswith("_1") else "Fabp3"
        target = -(CCONST[gene] / (KCAT_DLK[gene] * mult))

        reaction = m.reactions.get_by_id(rid)
        current = reaction.get_coefficient(pid)

        reaction.add_metabolites(
            {m.metabolites.get_by_id(pid): target - current},
            combine=True
        )

    # Block endogenous peroxidized-lipid production.
    for rid in HMR:
        m.reactions.get_by_id(rid).bounds = (0.0, 0.0)

    # Replace the presynaptic sink with a fixed C04717 source.
    m.reactions.get_by_id(HPODE_SOURCE).bounds = (0.0, 0.0)

    src = Reaction("SRC_PRE_LOAD")
    src.bounds = (load, load)
    src.add_metabolites({m.metabolites.get_by_id("PRE_C04717[c]"): 1.0})
    m.add_reactions([src])

    m.reactions.get_by_id(LINOLEATE_SINK).bounds = (0.0, 1000.0)

    # Apply the same diffusion-limited passive transport bound in WT and KO.
    cap = PASSIVE_F * FABP7_CAPACITY

    for rid in PASSIVE:
        if rid in have:
            m.reactions.get_by_id(rid).bounds = (0.0, cap)

    for rid in PASSIVE_BI:
        if rid in have:
            m.reactions.get_by_id(rid).bounds = (-cap, cap)

    # Apply abundance-derived protein-usage caps to all three FABPs.
    for gene, ratio in CAP_RATIO.items():
        if USAGE[gene] in have:
            m.reactions.get_by_id(USAGE[gene]).lower_bound = -(CAP_F7 * ratio)


rows = []


# Wild-type reference: fixed at the baseline DLKcat turnover numbers.
print()
print("Wild-type reference (FABP5/FABP3 at 1x DLKcat, FABP7 intact):")
print(f"{'load':>10} {'WT clearance':>14}")
print("-" * 26)

WT_REF = {}

for load in LOADS:
    with base as m:
        configure(m, WT_REF_MULT, load)
        clearance, _ = solve_max(m, DETOX)

    WT_REF[load] = clearance
    rows.append({
        "arm": "WT",
        "kcat_mult": WT_REF_MULT,
        "load": load,
        "clearance": clearance
    })

    print(f"{load:>10.5f} {clearance:>14.6f}", flush=True)


# Knockout: vary FABP5/FABP3 turnover while keeping WT fixed at 1x.
print()
print("Knockout, FABP5/FABP3 turnover multiplied (wild type held at 1x):")
print(f"{'mult':>7} {'load':>10} {'KO':>13} {'WT ref':>13} {'KO/WT':>8}")
print("-" * 55)

for mult in MULTS:
    for load in LOADS:
        with base as m:
            configure(m, mult, load)
            knock_out_model_genes(m, [FABP7_GENE])
            clearance, _ = solve_max(m, DETOX)

        rows.append({
            "arm": "KO",
            "kcat_mult": mult,
            "load": load,
            "clearance": clearance
        })

        wt = WT_REF[load]
        pct = 100.0 * clearance / wt if wt else float("nan")

        print(f"{mult:>7.4g} {load:>10.5f} {clearance:>13.6f} "
              f"{wt:>13.6f} {pct:>7.2f}%", flush=True)


df = pd.DataFrame(rows)


# Compare KO clearance with the fixed WT reference at each load.
wt_ref = (
    df[df.arm == "WT"][["load", "clearance"]]
    .rename(columns={"clearance": "WT_reference"})
)

wide = (
    df[df.arm == "KO"][["kcat_mult", "load", "clearance"]]
    .rename(columns={"clearance": "KO"})
    .merge(wt_ref, on="load", how="left")
)

wide["ko_over_wt_ref"] = wide["KO"] / wide["WT_reference"].replace(0, np.nan)
wide["pct_of_wild_type"] = wide["ko_over_wt_ref"] * 100


# Maximum clearance reached by each KO multiplier.
WT_PLATEAU = max(WT_REF.values())

plat = (
    df[df.arm == "KO"]
    .groupby("kcat_mult")["clearance"]
    .max()
    .rename("KO")
    .reset_index()
)

plat["WT_reference"] = WT_PLATEAU
plat["ko_over_wt_ref"] = plat["KO"] / WT_PLATEAU
plat["pct_of_wild_type"] = plat["ko_over_wt_ref"] * 100

ko1 = float(plat.loc[plat.kcat_mult == 1.0, "KO"].iloc[0])
plat["ko_rel_to_1x"] = plat["KO"] / ko1
FOLD = WT_PLATEAU / ko1


# Published column and sheet names.
COLUMNS = {
    "arm": "Genotype",
    "kcat_mult": "kcat multiplier",
    "load": "Delivered load",
    "clearance": "Clearance flux",
    "KO": "KO clearance flux",
    "WT_reference": "WT clearance flux",
    "pct_of_wild_type": "KO % of WT",
}

SHEET_COLUMNS = {
    "clearance_by_load": [
        "Genotype", "kcat multiplier", "Delivered load", "Clearance flux"
    ],
    "kcat_sensitivity_by_load": [
        "kcat multiplier", "Delivered load", "KO clearance flux",
        "WT clearance flux", "KO % of WT"
    ],
    "kcat_sensitivity_summary": [
        "kcat multiplier", "KO clearance flux",
        "WT clearance flux", "KO % of WT"
    ],
}


# Parameters used to calculate FABP transport capacities.
KCAT_ALL = dict(KCAT_DLK, Fabp7=6.2826)
CARRIER = {"Fabp7": "FABP7", "Fabp5": "FABP5", "Fabp3": "FABP3"}
UNIPROT = {"Fabp7": "E9Q0H6", "Fabp5": "Q3TLH6", "Fabp3": "Q5EBJ0"}

ENZYME_COEF = {
    "Fabp7": 0.911512,
    "Fabp5": 15.43627600192752,
    "Fabp3": 1.076659639709152
}

fabp_parameters = pd.DataFrame([
    {
        "Carrier": CARRIER[g],
        "UniProt": UNIPROT[g],
        "kcat DLKcat per s": KCAT_ALL[g],
        "Abundance relative to FABP7": CAP_RATIO[g],
        "Protein usage cap": round(CAP_F7 * CAP_RATIO[g], 9),
        "Enzyme coefficient": ENZYME_COEF[g],
        "Maximum C04717 transport capacity":
            CAP_F7 * CAP_RATIO[g] / ENZYME_COEF[g]
    }
    for g in ("Fabp7", "Fabp5", "Fabp3")
])


# Notes below the FABP parameter table.
FABP_NOTES = [
    "Notes",
    "Protein usage cap is a bound on the enzyme usage reaction "
    "(mmol·gDW⁻¹·h⁻¹).",
    "Maximum C04717 transport capacity = protein-usage cap "
    "÷ enzyme coefficient.",
]


def published(frame, sheet):
    return frame.rename(columns=COLUMNS)[SHEET_COLUMNS[sheet]]


def write_book(path):
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        fabp_parameters.to_excel(xl, sheet_name="fabp_parameters", index=False)

        for sheet, frame in (
            ("clearance_by_load", df),
            ("kcat_sensitivity_by_load", wide),
            ("kcat_sensitivity_summary", plat)
        ):
            published(frame, sheet).to_excel(xl, sheet_name=sheet, index=False)

        sheet = xl.book["fabp_parameters"]

        for cell in sheet[1]:
            cell.border = Border()

        first = len(fabp_parameters) + 3

        for i, line in enumerate(FABP_NOTES):
            cell = sheet.cell(row=first + i, column=1, value=line)

            if i == 0:
                cell.font = Font(bold=True)


# Save results.
RESULTS.mkdir(parents=True, exist_ok=True)
write_book(OUTPUT_XLSX)

print()
print(f"Saved: {OUTPUT_XLSX}")
print()
print("Compensation bound at plateau (wild type fixed at 1x):")
print(plat[["kcat_mult", "KO", "WT_reference", "pct_of_wild_type"]].to_string(index=False))
print()
print(f"Full compensation would need {FOLD:.1f}-fold higher FABP5/FABP3")
print(f"turnover than predicted ({math.log10(FOLD):.2f} log10 units).")