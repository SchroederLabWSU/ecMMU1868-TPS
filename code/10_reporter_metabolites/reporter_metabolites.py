"""
Reporter Metabolite Analysis — FABP7 KO vs WT (ecMMU1868-TPS)

Computes Patil-Nielsen reporter metabolite z-scores from WT and FABP7-KO
flux distributions under peroxidized-lipid stress. Reporter scores are
calculated from reaction-level flux changes, filtered to lipid-associated
metabolites.

Inputs:
  WT and FABP7-KO flux distributions from step 09
  model_S_triplets.csv
  model_rxn_meta.csv
  model_met_meta.csv

Output:
  reporter_metabolites_lipid.xlsx
"""


from pathlib import Path
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, Border, Side
from openpyxl.utils import get_column_letter

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
FABP7_KO_RESULTS = REPO_ROOT / "results" / "09_WT_vs_FABP7_ko_lipid_redox4"
REPORTER_RESULTS = REPO_ROOT / "results" / "10_reporter_metabolites"
OUTPUT_XLSX = REPORTER_RESULTS / "reporter_metabolites_lipid.xlsx"

WT = FABP7_KO_RESULTS  / "full_flux_wt_lipidperox.csv"
KO = FABP7_KO_RESULTS  / "full_flux_ko_lipidperox.csv"
S_TRIPLETS = REPORTER_RESULTS  / "model_S_triplets.csv"
RXN_META = REPORTER_RESULTS  / "model_rxn_meta.csv"
MET_META = REPORTER_RESULTS  / "model_met_meta.csv"



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
LIPID_SUBSYSTEMS = (
    "fatty acid", "lipid", "sphingo", "glycerolipid", "glycerophospholipid",
    "phosphatidylinositol", "cholesterol", "squalene", "steroid",
    "bile acid", "triacylglycerol", "carnitine shuttle",
    "linoleate", "arachidonic", "eicosanoid", "leukotriene",
    "ros detoxification", "glutathione",
)
COFACTOR_NAMES = {
    "proton", "water", "h2o", "oxygen", "o2", "co2", "carbon dioxide",
    "atp", "adp", "amp", "gtp", "gdp", "gmp",
    "nad", "nadh", "nadp", "nadph", "nadphx", "fad", "fadh2",
    "phosphate", "diphosphate", "pyrophosphate", "pi", "ppi", "coa",
    "ammonia", "ammonium", "nh3", "nh4", "hydrogen peroxide", "superoxide", "h2o2",
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
    "malonyl-coa", "malonylcoa", "acetoacetyl-coa", "hydroxymethylglutaryl-coa", "hmg-coa",
    "alanine", "aspartate", "glutamate", "glutamine", "glycine", "serine",
    "threonine", "valine", "leucine", "isoleucine", "acetyl-l-carnitine",
    "ribulose", "ribose", "xylulose", "sedoheptulose", "erythrose",
    "glycerol", "glycerol-3-phosphate", "glycerol 3-phosphate", "dihydroxyacetone",
    "maltotriose", "maltotetraose", "maltose", "starch", "glycogen",
}
GLYCOLYTIC_SUBSYSTEMS = (
    "glycolysis", "gluconeogenesis", "citric acid cycle", "tca cycle",
    "pyruvate metabolism", "pentose phosphate", "oxidative phosphorylation",
    "glyoxylate and dicarboxylate", "butanoate metabolism", "propanoate metabolism",
    "c5-branched dibasic acid", "urea cycle", "alanine and aspartate",
    "arginine and proline", "beta-alanine", "d-alanine", "glutamate metabolism",
    "glycine, serine, alanine, and threonine", "histidine metabolism",
    "lysine metabolism", "methionine and cysteine", "phenylalanine metabolism",
    "tryptophan metabolism", "tyrosine metabolism", "valine, leucine, and isoleucine",
    "taurine and hypotaurine", "nucleotide", "purine", "pyrimidine",
    "starch and sucrose", "fructose and mannose", "galactose", "aminosugar",
    "coa synthesis", "coa catabolism", "lipoate metabolism",
)
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

MAX_METABOLITE_DEGREE = 100 # maximum number of reactions a metabolite can be connected to and still be included in the analysis

CELL_TYPE_MAP = {"AST": "Astrocyte", "PRE": "Presynaptic neuron", "POST": "Postsynaptic neuron"}
COMPARTMENT_MAP = {
    "c": "cytosol", "e": "extracellular", "m": "mitochondrion",
    "r": "endoplasmic reticulum", "l": "lysosome", "g": "Golgi",
    "n": "nucleus", "x": "peroxisome", "syn": "synaptic cleft",
    "prec": "presynaptic cleft", "postc": "postsynaptic cleft",
    "pap": "perisynaptic astrocyte process", "u": "exchange pool", "q": "synaptic pool",
}


def strip_prefix_and_compartment(met_id):
    # Strip species prefix (AST_/PRE_/POST_) and compartment suffix, return lowercase base name.
    base = met_id
    for species_tag in ("AST_", "PRE_", "POST_"):
        if base.startswith(species_tag):
            base = base[len(species_tag):]
            break
    return base.split("[", 1)[0].lower() if "[" in base else base.lower()


def is_currency(met_id):
    # Return True if the metabolite is a common currency cofactor (ATP, NAD, water, etc.).
    return strip_prefix_and_compartment(met_id) in CURRENCY_IDS


def is_model_artifact(met_id):
    # Return True if the metabolite is a GECKO protein-pool or biomass pseudometabolite.
    base = strip_prefix_and_compartment(met_id)
    return (base.startswith("prot_") or base == "prot_pool"
            or base.startswith("biomass") or base.endswith("_biomass"))


def is_glycolytic(met_name, subsystem):
    # Return True if the metabolite belongs to central carbon or amino-acid metabolism.
    name = (met_name or "").lower().strip()
    subs = (subsystem or "").lower()
    return name in GLYCOLYTIC_NAMES or any(pathway in subs for pathway in GLYCOLYTIC_SUBSYSTEMS)


def is_lipid(met_name, subsystem):
    # Return True only for genuine lipid/fatty-acid metabolites.
    # Cofactors, glycolytic intermediates, and amino acids are excluded first.
    name = (met_name or "").lower().strip()
    subs = (subsystem or "").lower()
    if name in COFACTOR_NAMES:
        return False
    if is_glycolytic(name, subs):
        return False
    return any(keyword in name for keyword in LIPID_KEYWORDS) or \
           any(pathway in subs for pathway in LIPID_SUBSYSTEMS)


def parse_cell_type_and_compartment(met_id):
    # Extract cell type and compartment from a community metabolite ID.
    # e.g. AST_hdca[c] → ('Astrocyte', 'cytosol')
    cell_type = "Shared"
    body = met_id
    for species_tag in ("AST_", "PRE_", "POST_"):
        if met_id.startswith(species_tag):
            cell_type = CELL_TYPE_MAP[species_tag.rstrip("_")]
            body = met_id[len(species_tag):]
            break
    compartment = ""
    if "[" in body and "]" in body:
        comp_code = body.split("[")[-1].rstrip("]")
        compartment = COMPARTMENT_MAP.get(comp_code, comp_code)
    return cell_type, compartment


# Patil-Nielsen style reporter metabolite analysis
# For each reaction: compute |WT_flux - KO_flux| (magnitude of change) and
# (|WT_flux| - |KO_flux|) (signed direction: positive = suppressed in KO).
# Both are z-scored across all reactions so each metabolite's connected
# reactions contribute comparably regardless of absolute flux scale.

print("Loading flux tables .....", flush=True)
wt_fluxes = pd.read_csv(WT).rename(columns={"flux": "flux_wt"})
ko_fluxes = pd.read_csv(KO).rename(columns={"flux": "flux_ko"})
flux_table = wt_fluxes.merge(ko_fluxes, on="rxn_id", how="inner")

# Reaction-level flux change metrics for reporter scoring
flux_table["abs_change"]    = (flux_table["flux_wt"] - flux_table["flux_ko"]).abs()
flux_table["signed_change"] = flux_table["flux_wt"].abs() - flux_table["flux_ko"].abs()

# Standardise abs_change to z-scores across all reactions
delta_mean = flux_table["abs_change"].mean()
delta_std  = flux_table["abs_change"].std(ddof=0) or 1.0
flux_table["rxn_z_score"] = (flux_table["abs_change"] - delta_mean) / delta_std

# # Add reaction-level scores to each metabolite-reaction pair.
stoich = pd.read_csv(S_TRIPLETS)
stoich = stoich.join(
    flux_table.set_index("rxn_id")[["rxn_z_score", "signed_change"]],
    on="rxn_id"
).dropna(subset=["rxn_z_score"])

# Pre-filter metabolites by removing currency cofactors, GECKO artifacts, and
# highly connected hub metabolites (> MAX_METABOLITE_DEGREE connected reactions).
n_before_filter = stoich["met_id"].nunique()
stoich = stoich[~stoich["met_id"].map(is_currency)]
stoich = stoich[~stoich["met_id"].map(is_model_artifact)]
rxn_degree = stoich.groupby("met_id")["rxn_id"].nunique()
stoich = stoich[stoich["met_id"].isin(rxn_degree[rxn_degree <= MAX_METABOLITE_DEGREE].index)]
print(f"  metabolites: {n_before_filter} -> {stoich['met_id'].nunique()} "
      f"(after currency / artifact / degree filter)", flush=True)

# Aggregate reaction z-scores per metabolite: mean(z_j) * sqrt(n_j)
# Metabolites connected to more reactions get a higher score.
per_met_z  = stoich.groupby("met_id")["rxn_z_score"]
reporter_raw_score  = per_met_z.mean() * np.sqrt(per_met_z.count()); reporter_raw_score.name = "reporter_raw_score"
connected_reactions= per_met_z.count(); connected_reactions.name = "n_rxns"

# Convert raw reporter scores to background-normalized z-scores.
bg_mean = reporter_raw_score.mean()
bg_std = reporter_raw_score.std(ddof=0) or 1.0
reporter_z = (reporter_raw_score - bg_mean) / bg_std; reporter_z.name = "reporter_z"

# Signed direction: sum of signed_change per metabolite.
# Positive = metabolite's reactions carry more flux in WT than KO (suppressed in KO).
# Negative = metabolite's reactions carry more flux in KO than WT (induced in KO).
direction_score = stoich.groupby("met_id")["signed_change"].sum()
direction_score.name = "signed_delta"

# Annotate with metabolite names/formulas and subsystem labels
rxn_meta = pd.read_csv(RXN_META)
met_meta = pd.read_csv(MET_META)
stoich_with_subs = stoich.merge(rxn_meta[["rxn_id", "subsystem"]], on="rxn_id", how="left")
met_subsystems = (stoich_with_subs.groupby("met_id")["subsystem"]
        .apply(lambda v: ";".join(sorted({x for x in v.dropna().astype(str) if x}))))

results = (pd.concat([reporter_raw_score, connected_reactions, reporter_z, direction_score], axis=1)
        .join(met_meta.set_index("met_id"), how="left")
        .assign(subsystems=met_subsystems)
        .reset_index())

results["lipid_flag"] = [is_lipid(met_name, subsystem)
         for met_name, subsystem in zip(results["met_name"].fillna(""),
         results["subsystems"].fillna(""))]
results["direction"] = np.where(results["signed_delta"] > 0, "Suppressed",
         np.where(results["signed_delta"] < 0, "Induced", "Unchanged"))
results["abs_reporter_z"] = results["reporter_z"].abs()
results.sort_values("abs_reporter_z", ascending=False, inplace=True)

# Filter to lipid metabolites only
lipid_results = results[results["lipid_flag"]].copy()
lipid_results[["Cell Type", "Compartment"]] = lipid_results["met_id"].apply(lambda m: pd.Series(parse_cell_type_and_compartment(m)))
lipid_results = lipid_results.rename(columns={"met_name": "Metabolite Name", "n_rxns": "Connected Reactions",
         "direction": "Direction", "reporter_z": "Reporter score"})
lipid_results["Direction"] = lipid_results["Direction"].str.capitalize()

print(f"  lipid metabolites: {len(lipid_results)}", flush=True)

# Per-compartment view
full_per_comp = (lipid_results.sort_values("abs_reporter_z", ascending=False).reset_index(drop=True)
         .assign(**{"#": lambda d: range(1, len(d)+1)})
        [["#","Metabolite Name","Cell Type","Compartment","Reporter score","Connected Reactions","Direction"]])

# Collapse compartments, retaining the highest absolute reporter score for each metabolite within each cell type.
per_cell_type = (lipid_results.sort_values("abs_reporter_z", ascending=False)
         .drop_duplicates(subset=["Metabolite Name","Cell Type"], keep="first")
         .reset_index(drop=True))
per_cell_type["#"] = range(1, len(per_cell_type)+1)
per_cell_type = per_cell_type[["#","Metabolite Name","Cell Type","Compartment","Reporter score","Connected Reactions","Direction"]]

# Top 5% most extreme reporter scores (per cell type), directional only (Suppressed or Induced).
n_top = max(1, int(round(0.05 * len(per_cell_type))))
top5 = per_cell_type.head(n_top).copy()
top5_dir= top5[top5["Direction"].isin(["Suppressed","Induced"])].reset_index(drop=True)
top5_dir["#"] = range(1, len(top5_dir)+1)

print(f"  top 5% (per cell type): {len(top5)}  directional: {len(top5_dir)}", flush=True)
print("  Direction breakdown:", top5["Direction"].value_counts().to_dict(), flush=True)


# Excel export
HEADER_FONT = Font(name="Calibri", size=12, bold=True)
BODY_FONT = Font(name="Calibri", size=11)
THIN = Side(style="thin", color="BFBFBF")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
CENTER = Alignment(horizontal="center", vertical="center", wrap_text=True)
LEFT = Alignment(horizontal="left", vertical="center", wrap_text=True)
COL_WIDTHS = {"#": 5, "Metabolite Name": 44, "Cell Type": 20, "Compartment": 24,
               "Reporter score": 16, "Connected Reactions": 18, "Direction": 13}


def add_sheet(wb, name, df):
    ws = wb.create_sheet(name)
    cols = list(df.columns)
    for j, h in enumerate(cols, 1):
        c = ws.cell(row=1, column=j, value=h)
        c.font = HEADER_FONT; c.alignment = CENTER; c.border = BORDER
    for i, row in enumerate(df.itertuples(index=False), 2):
        for j, val in enumerate(row, 1):
            c = ws.cell(row=i, column=j, value=val)
            c.font = BODY_FONT; c.border = BORDER
            h = cols[j-1]
            if h == "Reporter score": c.number_format = "+0.000;-0.000;0.000"; c.alignment = CENTER
            elif h in ("#","Connected Reactions","Cell Type","Compartment","Direction"): c.alignment = CENTER
            else: c.alignment = LEFT
    for j, h in enumerate(cols, 1):
        ws.column_dimensions[get_column_letter(j)].width = COL_WIDTHS.get(h, 14)
    ws.row_dimensions[1].height = 28


def add_notes(wb):
    ws = wb.create_sheet("Notes", 0)
    ws.column_dimensions["A"].width = 90
    lines = [
        ("Lipid Reporter Metabolites — FABP7 KO vs WT (ecMMU1868-TPS)", True),
        ("", False),
        ("Analysis: Patil-Nielsen reporter z-scores from WT and FABP7-KO flux distributions",False),
        ("under forced peroxidized-linoleate input (C04717, 0.3 mmol/gDW/h).",False),
        ("Lipid-associated metabolites are identified by name and subsystem keywords.", False),
        ("Glycolytic, TCA, amino-acid, nucleotide, and short-chain-acid intermediates are excluded.", False),
        ("", False),
        ("Direction: Suppressed = lower summed reaction flux in FABP7 KO than WT.", False),
        ("           Induced = higher summed reaction flux in FABP7 KO than WT.", False),
        ("", False),
        ("Sheets:", True),
        ("  Top5pct_directional - Top 5% lipid rows (suppressed or induced only).", False),
        ("  Top5pct_all - Top 5% lipid rows (all directions).", False),
        ("  All_per_cell_type - All lipid rows, compartment-collapsed to highest |reporter score| per cell type.", False),
        ("  All_per_compartment - All lipid rows, one row per cell type × compartment (no collapsing).", False),
    ]
    for i, (txt, bold) in enumerate(lines, 1):
        c = ws.cell(row=i, column=1, value=txt)
        c.font = Font(name="Calibri", size=11, bold=bold)
        c.alignment = Alignment(wrap_text=True, vertical="top")


REPORTER_RESULTS.mkdir(parents=True, exist_ok=True)
wb = Workbook(); wb.remove(wb.active)
add_notes(wb)
add_sheet(wb, "Top5pct_directional", top5_dir)
add_sheet(wb, "Top5pct_all", top5)
add_sheet(wb, "All_per_cell_type", per_cell_type)
add_sheet(wb, "All_per_compartment", full_per_comp)
wb.save(OUTPUT_XLSX)
print(f"\nSaved: {OUTPUT_XLSX}")
