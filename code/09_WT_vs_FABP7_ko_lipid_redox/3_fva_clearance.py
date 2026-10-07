"""
FVA of the FABP7 Clearance Phenotype Across Alternative Optima

This script tests whether the WT and FABP7-KO clearance phenotype persists
across all alternative optimal flux distributions.

Clearance is measured at AST_LNLCOH_REPAIRc, whose flux represents total
clearance downstream of the five parallel glutathione-dependent reduction steps.

Input:  models/06_community_TPS/enzyme_constrained_tps/ecMMU1868-TPS.mat
Output: results/09_WT_vs_FABP7_ko_lipid_redox/fva_clearance.xlsx
"""

import sys
import time
import pandas as pd
from pathlib import Path
from cobra import Reaction
from cobra.flux_analysis import flux_variability_analysis

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
RESULTS = REPO_ROOT / "results" / "09_WT_vs_FABP7_ko_lipid_redox"
MAT_MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
OUTPUT_XLSX = RESULTS / "fva_clearance.xlsx"

# Load MATLAB v7.3 community model
sys.path.insert(0, str(SCRIPT_DIR))
from mat73_to_cobra import load_mat73_cobra


TOLERANCE = 1e-7
FORCE_FLUX = 1e-3
EPS_DETOX = 1e-4
FABP7_GENE = "12140"

# Clearance and C04717 reactions
SRC_ID = "SRC_PRE_C04717_LOAD"
PREC_MET = "PRE_C04717[prec]"
LINOLEATE_SINK = "AST_SK_lnlc_c"
CLEARANCE = "AST_LNLCOH_REPAIRc"
BIOMASS_AST = "AST_BIOMASS_reaction"
SELF_SYNTH = ["AST_HMR_2440_EXP_1", "AST_HMR_2440_EXP_2"]

# Passive C04717 transport routes
PROT_STD = "AST_prot_standard"
TD = [
    "AST_C04717td", "AST_C04717td_REV",
    "PRE_C04717td", "PRE_C04717td_REV",
    "POST_C04717td", "POST_C04717td_REV"
]

# Substrate-specific FABP5/FABP3 turnover parameters
PROT = {"Fabp5": "AST_prot_Q3TLH6", "Fabp3": "AST_prot_Q5EBJ0"}
KCAT = {"Fabp5": 0.2594, "Fabp3": 3.8233}
CCONST = {"Fabp5": 0.433019 * 9.2471, "Fabp3": 0.445155 * 9.2471}

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

# FABP protein-usage caps
USAGE_CAP = {
    "AST_usage_prot_E9Q0H6": 0.005823,
    "AST_usage_prot_Q3TLH6": 0.003140,
    "AST_usage_prot_Q5EBJ0": 0.000043
}


def set_cost(m, have, rxn_id, prot_id, new_cost):
    # Replace the enzyme coefficient for one transport reaction.
    if rxn_id not in have:
        return

    reaction = m.reactions.get_by_id(rxn_id)

    if prot_id not in {x.id for x in reaction.metabolites}:
        return

    met = m.metabolites.get_by_id(prot_id)
    reaction.add_metabolites(
        {met: -new_cost - reaction.get_coefficient(prot_id)},
        combine=True
    )


def build_model():
    # Build the corrected WT model used for the clearance analysis.
    print(f"loading {MAT_MODEL.name} ...", flush=True)

    m = load_mat73_cobra(str(MAT_MODEL))

    for solver in ("gurobi", "glpk"):
        try:
            m.solver = solver
            break
        except Exception:
            continue

    m.tolerance = TOLERANCE
    have = {r.id for r in m.reactions}
    mets = {x.id for x in m.metabolites}

    # Remove the prot_standard artifact from passive routes.
    if PROT_STD in mets:
        ps = m.metabolites.get_by_id(PROT_STD)

        for rid in TD:
            if rid in have:
                reaction = m.reactions.get_by_id(rid)

                if PROT_STD in {x.id for x in reaction.metabolites}:
                    reaction.add_metabolites(
                        {ps: -reaction.get_coefficient(PROT_STD)},
                        combine=True
                    )

    # Apply substrate-specific DLKcat costs.
    for rid in F5_RXNS:
        set_cost(m, have, rid, PROT["Fabp5"], CCONST["Fabp5"] / KCAT["Fabp5"])

    for rid in F3_RXNS:
        set_cost(m, have, rid, PROT["Fabp3"], CCONST["Fabp3"] / KCAT["Fabp3"])

    for rid in NEURON_NEW:
        if rid not in have:
            continue

        cell = rid.split("_")[0]
        gene = "Fabp5" if rid.endswith("_1") else "Fabp3"
        pid = f"{cell}_prot_{PROT[gene].rsplit('_', 1)[-1]}"

        if pid in mets:
            set_cost(m, have, rid, pid, CCONST[gene] / KCAT[gene])

    # Keep FABP5/FABP3 PAP-to-cytosol transport open.
    for rid in PAP_RXNS:
        if rid in have:
            m.reactions.get_by_id(rid).bounds = (0.0, 1000.0)

    # Add a fixed presynaptic 13-HPODE source.
    src = Reaction(SRC_ID)
    src.name = "Forced pre-synaptic 13-HPODE load"
    src.bounds = (FORCE_FLUX, FORCE_FLUX)
    src.add_metabolites({m.metabolites.get_by_id(PREC_MET): 1.0})
    m.add_reactions([src])

    if "PRE_SK_C04717_prec" in have:
        m.reactions.get_by_id("PRE_SK_C04717_prec").bounds = (0.0, 0.0)

    # Block endogenous astrocytic 13-HPODE production in both arms.
    for rid in SELF_SYNTH:
        if rid in have:
            m.reactions.get_by_id(rid).bounds = (0.0, 0.0)

    # Apply FABP protein-usage caps.
    for uid, cap in USAGE_CAP.items():
        if uid in have:
            m.reactions.get_by_id(uid).bounds = (-cap, 0.0)

    if LINOLEATE_SINK in have:
        m.reactions.get_by_id(LINOLEATE_SINK).upper_bound = 1000.0

    # Maximize biomass with a small secondary weight on clearance.
    objective = {m.reactions.get_by_id(BIOMASS_AST): 1.0}

    if CLEARANCE in have:
        objective[m.reactions.get_by_id(CLEARANCE)] = EPS_DETOX

    m.objective = objective

    print(f"  {len(m.reactions)} reactions, tolerance {m.tolerance:g}", flush=True)
    return m


def knockout(m):
    # Close reactions containing the FABP7 gene rule.
    blocked = [
        r.id for r in m.reactions
        if r.gene_reaction_rule and FABP7_GENE in r.gene_reaction_rule
    ]

    for rid in blocked:
        m.reactions.get_by_id(rid).bounds = (0.0, 0.0)

    return len(blocked)


def main():
    t0 = time.time()
    model = build_model()
    rows = []

    # Evaluate WT first, then apply the FABP7 knockout to the same model.
    for arm in ("WT", "KO"):
        if arm == "KO":
            print(f"\nKO: closed {knockout(model)} FABP7-dependent reactions",
                  flush=True)

        print(f"{arm}: objective {model.slim_optimize():.8f}", flush=True)

        # One process avoids worker re-import of this script.
        fva = flux_variability_analysis(
            model,
            reaction_list=[CLEARANCE],
            fraction_of_optimum=1.0,
            processes=1
        )

        low = float(fva["minimum"].iloc[0])
        high = float(fva["maximum"].iloc[0])

        rows.append({
            "Genotype": arm,
            "Clearance reaction": CLEARANCE,
            "Minimum flux": low,
            "Maximum flux": high
        })

        print(f"{arm}: clearance {low:.6g} to {high:.6g}", flush=True)

    clearance = pd.DataFrame(rows)

    RESULTS.mkdir(parents=True, exist_ok=True)
    clearance.to_excel(
        OUTPUT_XLSX,
        sheet_name="fva_clearance",
        index=False
    )

    # Non-overlap requires the KO maximum to remain below the WT minimum.
    wt, ko = clearance.iloc[0], clearance.iloc[1]
    result = (
        "disjoint"
        if ko["Maximum flux"] < wt["Minimum flux"]
        else "overlap"
    )

    print(f"\nKO maximum {ko['Maximum flux']:.6g} vs "
          f"WT minimum {wt['Minimum flux']:.6g}: {result}")
    print(f"wrote {OUTPUT_XLSX}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()