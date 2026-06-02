"""
Model Consistency Checks

Compares base model:iMM1865 and curated model:iMMU1867

Checks:
1. Stoichiometric consistency
2. Mass-imbalanced internal reactions
3. Charge-imbalanced internal reactions

Output:
consistency_checks_iMM1865_vs_iMMU1867.xlsx
"""


import pandas as pd
from pathlib import Path
from cobra.io import read_sbml_model
from memote.support.consistency import check_stoichiometric_consistency


# All paths are relative to the repository root
SCRIPT_DIR = Path(__file__).resolve().parent
MODELS_DIR = SCRIPT_DIR.parents[1] / "models"

MODEL_PATHS = {
    "iMM1865":  str(MODELS_DIR / "01_starting_model_iMM1865" / "iMM1865.xml"),
    "iMMU1867": str(MODELS_DIR / "02_curated_iMMU1867" / "iMMU1867.xml"),}

OUTPUT_XLSX = str(SCRIPT_DIR.parents[1] / "results" / "01_model_curation" / "consistency_checks_iMM1865_vs_iMMU1867.xlsx")
BIOMASS_ID = "BIOMASS_reaction"
TOLERANCE = 1e-6


def check_mass_and_charge_balance(model):
    """
    Identify mass and charge-imbalanced internal reactions.

    Boundary reactions (exchanges, skinks, and demands) and biomass reaction are excluded because
    they are expected to be unbalanced.
    """

    mass_imbalanced = []
    charge_imbalanced = []

    for reaction in model.reactions:

        if reaction.boundary or reaction.id == BIOMASS_ID:
            continue

        imbalance = reaction.check_mass_balance()

        if not imbalance:
            continue

        charge = imbalance.pop("charge", 0)

        if abs(charge) > TOLERANCE:
            charge_imbalanced.append({
                "Reaction_ID": reaction.id,
                "Reaction_Name": reaction.name,
                "Charge_Imbalance": charge,
                "Reaction": reaction.reaction,
            })

        if imbalance:
            mass_imbalanced.append({
                "Reaction_ID": reaction.id,
                "Reaction_Name": reaction.name,
                "Mass_Imbalance": "; ".join(
                    f"{element}: {value:g}"
                    for element, value in imbalance.items()
                ),
                "Reaction": reaction.reaction,
            })

    return pd.DataFrame(mass_imbalanced), pd.DataFrame(charge_imbalanced)


def check_model_stoichiometry(model):
    # Check stoichiometric consistency using MEMOTE functions.

    is_consistent = check_stoichiometric_consistency(model)

    return is_consistent


def analyze_model(model_name, model_path):
    # Load model and run all consistency checks.

    print(f"\nAnalyzing {model_name}")
    print("-----" * 5)

    model = read_sbml_model(model_path)

    print(f"Reactions:   {len(model.reactions)}")
    print(f"Metabolites: {len(model.metabolites)}")
    print(f"Genes:       {len(model.genes)}")

    print("Checking mass and charge balance...")
    mass_df, charge_df = check_mass_and_charge_balance(model)

    print("Checking stoichiometric consistency...")
    stoich_consistent = check_model_stoichiometry(model)

    summary = {
        "Model": model_name,
        "Reactions": len(model.reactions),
        "Metabolites": len(model.metabolites),
        "Genes": len(model.genes),
        "Stoichiometrically consistent": "Yes" if stoich_consistent else "No",
        "Mass-imbalanced reactions": len(mass_df),
        "Charge-imbalanced reactions": len(charge_df),
    }

    result_tables = {
        "Mass_Imbalanced_Reactions": mass_df,
        "Charge_Imbalanced_Reactions": charge_df,
    }

    return summary, result_tables


# Main analysis 

all_summaries = []
all_results = {}

for model_name, model_path in MODEL_PATHS.items():

    summary, result_tables = analyze_model(model_name, model_path)
    all_summaries.append(summary)

    for table_name, table in result_tables.items():
        sheet_name = f"{model_name}_{table_name}"[:31]
        all_results[sheet_name] = table

summary_df = pd.DataFrame(all_summaries)


# Save results to Excel

print(f"\nWriting results to {OUTPUT_XLSX}")
Path(OUTPUT_XLSX).parent.mkdir(parents=True, exist_ok=True)

with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:

    summary_df.to_excel(
        writer,
        sheet_name="Summary",
        index=False
    )

    for sheet_name, table in all_results.items():

        if table.empty:
            table = pd.DataFrame({"Note": ["None found"]})

        table.to_excel(
            writer,
            sheet_name=sheet_name,
            index=False
        )

print("Done.")