"""
Reporter-Metabolite Analysis of FABP7 KO vs WT

Computes reporter-metabolite scores from WT and FABP7-KO flux distributions
under forced peroxidized-lipid stress across 100 turnover-number
parameterizations.

Each parameterization is solved under both ATP-maintenance and biomass
objectives using the same kcat draw.

Inputs:
  models/06_community_TPS/enzyme_constrained_tps/ecMMU1868-TPS.mat
  model_S_triplets.csv
  model_met_meta.csv
  kcat_parameters.xlsx

Output:
  reporter_metabolites_bootstrap.xlsx
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from cobra import Reaction
from cobra.flux_analysis import pfba

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
REPORTER_RESULTS = REPO_ROOT / "results" / "10_reporter_metabolites"
OUTPUT_XLSX = REPORTER_RESULTS / "reporter_metabolites_bootstrap.xlsx"

MODEL = REPO_ROOT / "models" / "06_community_TPS" / "enzyme_constrained_tps" / "ecMMU1868-TPS.mat"
S_TRIPLETS = REPORTER_RESULTS / "model_S_triplets.csv"
MET_META = REPORTER_RESULTS / "model_met_meta.csv"
KCAT_PARAMS = REPORTER_RESULTS / "kcat_parameters.xlsx"

# Load MATLAB v7.3 community model
sys.path.insert(0, str(SCRIPT_DIR))
from mat73_to_cobra import load_mat73_cobra

# Bootstrap settings
N_RUNS = 100
BASE_SEED = 20260818
P_DLKCAT = 0.5
OBJECTIVES = ("atpm", "biomass")

# Solver and reporter settings
TOLERANCE = 1e-7
GROWTH_FRACTION = 0.99
NEAR_ZERO = 1e-6
TOP_N = 10         # Threshold for counting entries within the top ten
MAX_DEGREE = 100   # Exclude metabolites connected to more than 100 reactions

# Forced peroxidized-lipid stress
FORCE_FLUX = 1e-3
FABP7_GENE = "12140"

SRC_ID = "SRC_PRE_C04717_LOAD"
PREC_MET = "PRE_C04717[prec]"
LINOLEATE_SINK = "AST_SK_lnlc_c"
LIPID_REPAIR = "AST_LNLCOH_REPAIRc"

BIOMASS_AST = "AST_BIOMASS_reaction"
BIOMASS_ALL = [
    "AST_BIOMASS_reaction",
    "PRE_BIOMASS_reaction",
    "POST_BIOMASS_reaction"
]
ATPM_ALL = ["AST_ATPM", "PRE_ATPM", "POST_ATPM"]

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


# Terms used to identify lipid-associated metabolites
LIPID_KEYWORDS = (
    "fatty acid", "fatty acyl", "lipid",
    "sphing", "ceramide", "phosphatidyl", "lysophosphatidyl", "lyso",
    "linole", "palmit", "stearoyl", "stearate", "oleoyl", "oleate",
    "arachid", "eicos", "docosa", "prostaglandin", "leukotriene",
    "cholester", "lanost", "sterol", "bile acid", "bile-acid",
    "triacylglycerol", "diacylglycerol", "monoacylglycerol",
    "carnitine", "acylcoa", "acyl-coa",
    "long-chain acyl", "medium-chain acyl", "short-chain acyl",
    "lnlc", "hdca", "ocdca", "stcoa", "pmtcoa", "c04717",
)

COFACTOR_NAMES = {
    "proton", "water", "h2o", "oxygen", "o2", "co2", "carbon dioxide",
    "atp", "adp", "amp", "gtp", "gdp", "gmp",
    "nad", "nadh", "nadp", "nadph", "nadphx", "fad", "fadh2",
    "phosphate", "diphosphate", "pyrophosphate", "pi", "ppi", "coa",
    "ammonia", "ammonium", "nh3", "nh4",
    "hydrogen peroxide", "superoxide", "h2o2",
}

GLYCOLYTIC_NAMES = {
    "glucose", "glucose-6-phosphate", "fructose-6-phosphate",
    "fructose-1,6-bisphosphate", "fructose 1,6-bisphosphate",
    "dihydroxyacetone phosphate", "dhap",
    "glyceraldehyde-3-phosphate", "glyceraldehyde 3-phosphate", "g3p",
    "1,3-bisphosphoglycerate", "3-phosphoglycerate", "2-phosphoglycerate",
    "phosphoenolpyruvate", "pep", "pyruvate", "lactate", "l-lactate", "d-lactate",
    "acetyl-coa", "acetyl coa", "acetylcoa", "acetyl coenzyme a",
    "oxaloacetate", "citrate", "isocitrate", "alpha-ketoglutarate", "2-oxoglutarate",
    "succinyl-coa", "succinylcoa", "succinate", "fumarate", "malate", "l-malate",
    "malonyl-coa", "malonylcoa", "acetoacetyl-coa",
    "hydroxymethylglutaryl-coa", "hmg-coa",
    "alanine", "aspartate", "glutamate", "glutamine", "glycine", "serine",
    "threonine", "valine", "leucine", "isoleucine", "acetyl-l-carnitine",
    "ribulose", "ribose", "xylulose", "sedoheptulose", "erythrose",
    "glycerol", "glycerol-3-phosphate", "glycerol 3-phosphate",
    "dihydroxyacetone", "maltotriose", "maltotetraose", "maltose",
    "starch", "glycogen",
}

CURRENCY_IDS = {
    "h", "h2o", "o2", "co2", "hco3", "o2s", "h2o2",
    "na1", "k", "cl", "ca2", "mg2", "fe2", "fe3", "cu", "cu2",
    "zn2", "mn2", "cobalt2", "mobd", "so4", "so3", "no", "no2",
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "utp", "udp", "ump",
    "ctp", "cdp", "cmp", "itp", "idp", "imp",
    "datp", "dadp", "damp", "dgtp", "dgdp", "dgmp",
    "dctp", "dcdp", "dcmp", "dttp", "dtdp", "dtmp",
    "nad", "nadh", "nadp", "nadph", "nadphx", "fad", "fadh2",
    "fmn", "fmnh2", "q", "qh2", "q10", "q10h2", "mqn8", "mql8",
    "trdox", "trdrd", "thrdx", "thrdxox",
    "pi", "ppi", "pppi", "coa", "nh3", "nh4",
}

CELL_TYPE_MAP = {
    "AST": "Astrocyte",
    "PRE": "Presynaptic neuron",
    "POST": "Postsynaptic neuron"
}

COMPARTMENT_MAP = {
    "c": "cytosol",
    "e": "extracellular",
    "m": "mitochondrion",
    "r": "endoplasmic reticulum",
    "l": "lysosome",
    "g": "Golgi",
    "n": "nucleus",
    "x": "peroxisome",
    "syn": "synaptic cleft",
    "prec": "presynaptic cleft",
    "postc": "postsynaptic cleft",
    "pap": "perisynaptic astrocyte process",
    "u": "exchange pool",
    "q": "synaptic pool",
}

# Cleaner names used only in the output.
NICE_NAMES = {"Linoleic Coenzyme A": "Linoleoyl-CoA"}
NICE_CELLS = {
    "Presynaptic neuron": "Pre-synaptic neuron",
    "Postsynaptic neuron": "Post-synaptic neuron"
}


def set_cost(m, have, rxn_id, prot_id, new_cost):
    # Replace the protein coefficient for one reaction.
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


def build_base():
    # Load the model and apply corrections shared by every parameterization.
    print(f"loading {MODEL.name} ...", flush=True)

    m = load_mat73_cobra(str(MODEL))

    for solver in ("gurobi", "glpk"):
        try:
            m.solver = solver
            break
        except Exception:
            continue

    m.tolerance = TOLERANCE
    have = {r.id for r in m.reactions}
    mets = {x.id for x in m.metabolites}

    print(f"  {len(m.reactions)} rxns | {len(m.metabolites)} mets", flush=True)

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
        set_cost(m, have, rid, PROT["Fabp5"],
                 CCONST["Fabp5"] / KCAT["Fabp5"])

    for rid in F3_RXNS:
        set_cost(m, have, rid, PROT["Fabp3"],
                 CCONST["Fabp3"] / KCAT["Fabp3"])

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
    assert PREC_MET in mets, f"{PREC_MET} not in model"

    src = Reaction(SRC_ID)
    src.name = "Forced pre-synaptic 13-HPODE load"
    src.bounds = (FORCE_FLUX, FORCE_FLUX)
    src.add_metabolites({m.metabolites.get_by_id(PREC_MET): 1.0})
    m.add_reactions([src])

    if "PRE_SK_C04717_prec" in have:
        m.reactions.get_by_id("PRE_SK_C04717_prec").bounds = (0.0, 0.0)

    # close astrocytic self-synthesis.
    for rid in SELF_SYNTH:
        if rid in have:
            m.reactions.get_by_id(rid).bounds = (0.0, 0.0)

    # Apply FABP protein-usage caps.
    for uid, cap in USAGE_CAP.items():
        if uid in have:
            m.reactions.get_by_id(uid).bounds = (-cap, 0.0)

    if LINOLEATE_SINK in have:
        m.reactions.get_by_id(LINOLEATE_SINK).upper_bound = 1000.0

    print(f"  base built: {len(m.reactions)} reactions, "
          f"tolerance {m.tolerance:g}", flush=True)

    return m


def biomass_mu_star(m):
    # Maximum WT growth for the current kcat parameterization.
    m.objective = {m.reactions.get_by_id(BIOMASS_AST): 1.0}
    m.objective_direction = "max"
    return float(m.slim_optimize())


def apply_objective(m, objective, growth_floor=None):
    # Set the physiological constraint and maximize astrocytic lipid repair.
    have = {r.id for r in m.reactions}

    if objective == "atpm":
        # Block growth and fix each ATPM reaction at its model lower bound.
        for rid in BIOMASS_ALL:
            if rid in have:
                m.reactions.get_by_id(rid).bounds = (0.0, 0.0)

        for rid in ATPM_ALL:
            if rid in have:
                reaction = m.reactions.get_by_id(rid)
                reaction.bounds = (reaction.lower_bound, reaction.lower_bound)

    else:
        # Maintain growth at at least 99% of the WT optimum for this run.
        m.reactions.get_by_id(BIOMASS_AST).bounds = (growth_floor, 1000.0)

    m.objective = {m.reactions.get_by_id(LIPID_REPAIR): 1.0}
    m.objective_direction = "max"


def knockout(m):
    # Close every reaction whose GPR contains FABP7.
    for reaction in [
        r for r in m.reactions
        if r.gene_reaction_rule and FABP7_GENE in r.gene_reaction_rule
    ]:
        reaction.bounds = (0.0, 0.0)


def solve_pfba(m, tag):
    # pFBA selects one parsimonious solution at the exact optimum.
    try:
        return pfba(m, fraction_of_optimum=1.0)
    except Exception as e:
        print(f"    {tag}: pFBA failed ({type(e).__name__})", flush=True)
        return None


def load_pairs():
    # Load reactions with both BRENDA and DLKcat turnover numbers.
    df = pd.read_excel(KCAT_PARAMS, sheet_name=0)
    df = df.rename(columns={
        c: c.strip().replace(" ", "_").lower()
        for c in df.columns
    })

    paired = df[
        (df.source_assigned.astype(str) == "brenda")
        & df.kcat_dlkcat_cached.notna()
        & (df.kcat_dlkcat_cached > 0)
        & (df.kcat_assigned > 0)
    ].copy()

    # Remove pairs that agree within rounding.
    paired = paired[
        (paired.kcat_assigned - paired.kcat_dlkcat_cached).abs()
        / paired.kcat_assigned > 1e-6
    ]

    return {
        str(r.reaction): (
            float(r.kcat_assigned),
            float(r.kcat_dlkcat_cached)
        )
        for r in paired.itertuples(index=False)
    }


def cache_coefficients(m, pairs):
    # Cache the baseline enzyme coefficient for each resampleable reaction.
    have = {r.id for r in m.reactions}
    coef = {}

    for rid in pairs:
        if rid not in have:
            continue

        reaction = m.reactions.get_by_id(rid)

        for met in reaction.metabolites:
            if "_prot_" in met.id:
                coef[rid] = (met.id, reaction.get_coefficient(met.id))
                break

    print(f"  paired in sheet {len(pairs)} | "
          f"resampleable in model {len(coef)}", flush=True)

    return coef


def apply_draw(m, coef, pairs, take_dlkcat):
    # Rescale coefficients only for reactions assigned DLKcat.
    for (rid, (pid, original)), use_dlkcat in zip(coef.items(), take_dlkcat):
        if not use_dlkcat:
            continue

        kcat_assigned, kcat_dlkcat = pairs[rid]
        target = original * (kcat_assigned / kcat_dlkcat)

        reaction = m.reactions.get_by_id(rid)
        shift = target - reaction.get_coefficient(pid)

        if abs(shift) > 1e-15:
            reaction.add_metabolites(
                {m.metabolites.get_by_id(pid): shift},
                combine=True
            )


def met_base(met_id):
    # Remove the cell prefix and compartment suffix.
    base = met_id

    for tag in ("AST_", "PRE_", "POST_"):
        if base.startswith(tag):
            base = base[len(tag):]
            break

    return base.split("[", 1)[0].lower() if "[" in base else base.lower()


def is_currency(met_id):
    return met_base(met_id) in CURRENCY_IDS


def is_model_artifact(met_id):
    # Exclude GECKO protein-pool and biomass pseudo-metabolites.
    base = met_base(met_id)

    return (
        base.startswith("prot_")
        or base == "prot_pool"
        or base.startswith("biomass")
        or base.endswith("_biomass")
    )


def is_lipid(met_name):
    # Keep lipid species while excluding cofactors and central-carbon metabolites.
    name = (met_name or "").lower().strip()

    if name in COFACTOR_NAMES or name in GLYCOLYTIC_NAMES:
        return False

    return any(keyword in name for keyword in LIPID_KEYWORDS)


def cell_and_compartment(met_id):
    # Parse cell type and compartment from the metabolite ID.
    cell, body = "Shared", met_id

    for tag in ("AST_", "PRE_", "POST_"):
        if met_id.startswith(tag):
            cell = CELL_TYPE_MAP[tag.rstrip("_")]
            body = met_id[len(tag):]
            break

    if "[" not in body or "]" not in body:
        return cell, ""

    code = body.split("[")[-1].rstrip("]")
    return cell, COMPARTMENT_MAP.get(code, code)


def safe_std(series):
    # Replace zero or invalid standard deviations with 1.
    std = float(series.std(ddof=0))
    return 1.0 if (not np.isfinite(std) or std == 0.0) else std


def score_reporters(flux_wt, flux_ko, triplets, met_names):
    # Calculate reporter scores from absolute WT-KO reaction-flux differences.
    flux = (
        flux_wt.rename(columns={"flux": "flux_wt"})
        .merge(flux_ko.rename(columns={"flux": "flux_ko"}), on="rxn_id")
    )

    flux["abs_change"] = (flux.flux_wt - flux.flux_ko).abs()
    flux["signed_change"] = flux.flux_wt.abs() - flux.flux_ko.abs()
    flux["rxn_z"] = (
        (flux.abs_change - flux.abs_change.mean())
        / safe_std(flux.abs_change)
    )

    stoich = (
        triplets
        .join(
            flux.set_index("rxn_id")[["rxn_z", "signed_change"]],
            on="rxn_id"
        )
        .dropna(subset=["rxn_z"])
    )

    # Remove currency metabolites, model artifacts, and highly connected hubs.
    stoich = stoich[~stoich.met_id.map(is_currency)]
    stoich = stoich[~stoich.met_id.map(is_model_artifact)]

    degree = stoich.groupby("met_id")["rxn_id"].nunique()
    stoich = stoich[
        stoich.met_id.isin(degree[degree <= MAX_DEGREE].index)
    ]

    per_met = stoich.groupby("met_id")["rxn_z"]
    raw = per_met.mean() * np.sqrt(per_met.count())

    # Positive signed_delta indicates lower summed absolute flux in KO.
    scores = pd.DataFrame({
        "met_id": raw.index,
        "reporter_z": ((raw - raw.mean()) / safe_std(raw)).values,
        "signed_delta": stoich.groupby("met_id")["signed_change"].sum().values,
        "n_rxns": per_met.count().values,
    })

    scores["name"] = scores.met_id.astype(str).map(met_names)

    parsed = scores.met_id.astype(str).map(cell_and_compartment)
    scores["cell"] = [p[0] for p in parsed]
    scores["compartment"] = [p[1] for p in parsed]

    return scores


def rank_one_run(scores):
    # Rank lipid metabolite-cell combinations for one parameterization.
    text = scores.name.astype(str)

    is_lip = pd.Series(
        [is_lipid(name) for name in scores.name.fillna("")],
        index=scores.index
    )

    # Exclude shared/cleft pools, bulk pseudo-metabolites, and the forced C04717 species.
    drop = (
        (scores.cell == "Shared")
        | scores.compartment.astype(str).str.contains(
            "cleft", case=False, na=False
        )
        | text.str.contains("bulk cytosolic pool", case=False, na=False)
        | text.str.contains("hydroperoxylinoleic", case=False, na=False)
    )

    lipids = scores[
        is_lip & scores.name.notna() & ~drop
    ].copy()

    lipids["entry"] = (
        lipids.name.astype(str) + " | " + lipids.cell.astype(str)
    )
    lipids["size"] = lipids.reporter_z.abs()

    # Keep the strongest compartment for each metabolite-cell combination.
    lipids = (
        lipids
        .sort_values("size", ascending=False)
        .drop_duplicates(subset=["entry"], keep="first")
    )

    # Near-zero WT-KO changes do not enter the ranking.
    ranked = lipids[
        lipids.signed_delta.abs() > NEAR_ZERO
    ].copy()

    ranked["rank"] = ranked["size"].rank(
        ascending=False,
        method="first"
    )

    return ranked.sort_values("rank")


def summarise(ranked):
    # Summarize ranking frequency and direction across all parameterizations.
    by_entry = ranked.groupby("entry")

    table = pd.DataFrame({
        "Median reporter score": by_entry["reporter_z"].median().abs(),
        "Ranked first (of 100)": by_entry["rank"].apply(
            lambda r: int((r == 1).sum())
        ),
        "In top three (of 100)": by_entry["rank"].apply(
            lambda r: int((r <= 3).sum())
        ),
        "Within top ten (of 100)": by_entry["rank"].apply(
            lambda r: int((r <= TOP_N).sum())
        ),
        "Parameterizations scored": by_entry["run"].nunique(),
        "Median rank when scored": by_entry["rank"].median(),
        "Times suppressed": by_entry["signed_delta"].apply(
            lambda s: int((s > 0).sum())
        ),
    }).reset_index()

    # Direction consistency is calculated among parameterizations where scored.
    scored = table["Parameterizations scored"]

    table["Direction"] = np.where(
        table["Times suppressed"] >= scored / 2,
        "Suppressed",
        "Induced"
    )

    table["Direction consistency"] = (
        np.maximum(
            table["Times suppressed"],
            scored - table["Times suppressed"]
        ) / scored
    ).round(3)

    table[["Metabolite", "Cell type"]] = table.entry.str.split(
        " | ",
        regex=False,
        expand=True
    )

    table["Metabolite"] = table["Metabolite"].replace(NICE_NAMES)
    table["Cell type"] = table["Cell type"].replace(NICE_CELLS)

    # Order by how often each entry ranked first, then by supporting ranks.
    table = table.sort_values(
        [
            "Ranked first (of 100)",
            "In top three (of 100)",
            "Within top ten (of 100)",
            "Median rank when scored"
        ],
        ascending=[False, False, False, True]
    ).reset_index(drop=True)

    table.insert(0, "Rank", range(1, len(table) + 1))

    return table[[
        "Rank", "Metabolite", "Cell type", "Median reporter score",
        "Ranked first (of 100)", "In top three (of 100)",
        "Within top ten (of 100)", "Parameterizations scored",
        "Median rank when scored", "Direction", "Direction consistency"
    ]]


def per_run_sheet(ranked):
    # Store every ranked entry so the summary can be reproduced.
    out = pd.DataFrame({
        "Parameterization": ranked.run.astype(int),
        "Rank": ranked["rank"].astype(int),
        "Metabolite": ranked.name.replace(NICE_NAMES),
        "Cell type": ranked.cell.replace(NICE_CELLS),
        "Compartment": ranked.compartment,
        "Reporter score": ranked.reporter_z,
        "WT - KO flux difference": ranked.signed_delta,
        "Connected reactions": ranked.n_rxns,
        "Direction": np.where(
            ranked.signed_delta > 0,
            "Suppressed",
            "Induced"
        ),
    })

    return out.sort_values(
        ["Parameterization", "Rank"]
    ).reset_index(drop=True)


# Notes included in the output workbook.
NOTES = pd.DataFrame([
    (
        "Overview",
        "Reporter-metabolite results comparing FABP7 knockout with wild type in "
        "ecMMU1868-TPS under forced peroxidized-lipid load (C04717, 0.001 "
        "mmol·gDW⁻¹·h⁻¹) across 100 kcat bootstrap parameterizations."
    ),
    (
        "Parameterizations",
        "Reactions with both BRENDA and DLKcat values were resampled across 100 "
        "parameterizations. The same parameterizations were used for the ATPM "
        "and biomass objectives."
    ),
    (
        "ATP-maintenance objective",
        "Growth was blocked, ATP-maintenance demand was imposed, and astrocytic "
        "clearance was maximized. This is the primary objective in the manuscript."
    ),
    (
        "Biomass objective",
        "Community growth was maintained at ≥99% of its optimum while "
        "astrocytic clearance was maximized."
    ),
    (
        "Reporter score",
        "Patil–Nielsen scores are standardized within each parameterization "
        "and indicate magnitude only; ranks are compared across parameterizations."
    ),
    (
        "Counts",
        "Counts are out of 100 parameterizations; scored indicates the number "
        "in which an entry entered the ranking."
    ),
    (
        "Direction",
        "Suppressed or induced indicates decreased or increased summed absolute "
        "flux in the knockout; consistency is the fraction of scored "
        "parameterizations with that direction."
    ),
    (
        "Exclusions",
        "Shared pools, cleft metabolites, and C04717 were excluded from reporter "
        "ranking."
    ),
], columns=["Item", "Detail"])


def write_workbook(summary, per_run):
    # Write notes, summaries, and per-parameterization results.
    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as book:
        NOTES.to_excel(book, sheet_name="Notes", index=False)

        summary["atpm"].to_excel(
            book, sheet_name="ATPM_summary", index=False
        )
        summary["biomass"].to_excel(
            book, sheet_name="biomass_summary", index=False
        )

        per_run["atpm"].to_excel(
            book,
            sheet_name="ATPM_per_parameterization",
            index=False
        )
        per_run["biomass"].to_excel(
            book,
            sheet_name="biomass_per_parameterization",
            index=False
        )

        # Set readable column widths and freeze the header row.
        for sheet in book.book.worksheets:
            for column in sheet.columns:
                longest = max(
                    len(str(cell.value or ""))
                    for cell in column
                )
                sheet.column_dimensions[
                    column[0].column_letter
                ].width = min(max(longest + 2, 12), 95)

            sheet.freeze_panes = "A2"


def main():
    REPORTER_RESULTS.mkdir(parents=True, exist_ok=True)

    triplets = pd.read_csv(S_TRIPLETS)
    met_names = (
        pd.read_csv(MET_META)
        .set_index("met_id")["met_name"]
    )

    print(f"metadata: {len(triplets)} triplets, "
          f"{len(met_names)} metabolites", flush=True)

    base = build_base()
    pairs = load_pairs()
    coef = cache_coefficients(base, pairs)

    collected = {objective: [] for objective in OBJECTIVES}

    for run in range(N_RUNS):
        # Use the same deterministic kcat draw for both objectives.
        seed = BASE_SEED + run
        take = (
            np.random.default_rng(seed).random(len(coef))
            < P_DLKCAT
        )

        with base as m:
            apply_draw(m, coef, pairs, take)

            print(f"run {run:3d}  seed {seed}  DLKcat draws "
                  f"{int(take.sum())}/{len(coef)}", flush=True)

            for objective in OBJECTIVES:
                growth_floor = None

                if objective == "biomass":
                    with m as mstar:
                        growth_floor = (
                            GROWTH_FRACTION
                            * biomass_mu_star(mstar)
                        )

                with m as mo:
                    apply_objective(
                        mo,
                        objective,
                        growth_floor
                    )

                    wt = solve_pfba(
                        mo,
                        f"{objective}/WT"
                    )

                    with mo as mk:
                        knockout(mk)
                        ko = solve_pfba(
                            mk,
                            f"{objective}/KO"
                        )

                if wt is None or ko is None:
                    print(f"  {objective}: dropped", flush=True)
                    continue

                scores = score_reporters(
                    wt.fluxes.rename_axis("rxn_id").reset_index(name="flux"),
                    ko.fluxes.rename_axis("rxn_id").reset_index(name="flux"),
                    triplets,
                    met_names
                )

                ranked = rank_one_run(scores)
                ranked["run"] = run

                collected[objective].append(
                    ranked[[
                        "run", "rank", "name", "cell", "compartment",
                        "reporter_z", "signed_delta", "n_rxns", "entry"
                    ]]
                )

                print(
                    f"  {objective}: {len(ranked)} ranked entries, "
                    f"top = {ranked.iloc[0]['name']} "
                    f"({ranked.iloc[0].cell})",
                    flush=True
                )

    summary, per_run = {}, {}

    for objective in OBJECTIVES:
        ranked = pd.concat(
            collected[objective],
            ignore_index=True
        )

        summary[objective] = summarise(ranked)
        per_run[objective] = per_run_sheet(ranked)

        print(f"{objective:8s} {ranked.run.nunique()} parameterizations | "
              f"{len(summary[objective])} entries", flush=True)

    write_workbook(summary, per_run)
    print(f"\nwrote {OUTPUT_XLSX}")


if __name__ == "__main__":
    main()