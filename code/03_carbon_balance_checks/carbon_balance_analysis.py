'''
Carbon Balance Analysis

Carbon balance analysis evaluates whether carbon entering
the metabolic network is properly accounted for through
metabolite secretion and biomass production.

Runs pFBA on both iMM1865 and iMMU1867 and exports results for each model
into separate sheets in the same Excel workbook. This sohws the
before-and-after curation comparison in Table 1.

This script:
1. Loads each SBML model
2. Runs parsimonious Flux Balance Analysis (pFBA) - minimizes total rxn flux while
maximizing biomass, giving a more realistic flux distribution for carbon balance analysis.
3. Calculates carbon uptake, secretion, and biomass carbon drain
4. Exports results to Excel

Output: carbon_balance_results.xlsx in results/03_carbon_balance_checks/
'''


import re
import cobra
import pandas as pd
from pathlib import Path
from cobra.flux_analysis import pfba


# All paths are relative to the repository root.

SCRIPT_DIR = Path(__file__).resolve().parent
MODELS_DIR = SCRIPT_DIR.parents[1] / "models"

MODELS = {
    "iMM1865":  MODELS_DIR / "01_starting_model_iMM1865" / "iMM1865.xml",
    "iMMU1867": MODELS_DIR / "02_curated_iMMU1867"       / "iMMU1867.xml",
}

OUTPUT_XLSX = SCRIPT_DIR.parents[1] / "results" / "03_carbon_balance_checks" / "carbon_balance_results.xlsx"

BIOMASS_ID = "BIOMASS_reaction"
EXTRACELLULAR_COMPARTMENT = "e"
FLUX_CUTOFF = 1e-6


def count_carbons(formula):
    # Return the number of carbon atoms in a metabolite formula.
    if not formula:
        return 0

    match = re.search(r"C(\d*)", formula)

    if not match:
        return 0

    carbon_number = match.group(1)

    if carbon_number == "":
        return 1

    return int(carbon_number)


def get_exchange_metabolite(reaction):
    # Return the single metabolite in an exchange reaction.
    metabolites = list(reaction.metabolites.keys())

    if len(metabolites) != 1:
        return None

    return metabolites[0]


def collect_exchange_carbon_fluxes(model, solution):
    # Calculate carbon uptake and secretion through exchange reactions.
    carbon_inputs = []
    carbon_outputs = []
    skipped = []

    for reaction in model.exchanges:
        flux = solution.fluxes.get(reaction.id, 0)

        if abs(flux) < FLUX_CUTOFF:
            continue

        metabolite = get_exchange_metabolite(reaction)

        if metabolite is None:
            continue

        if metabolite.compartment != EXTRACELLULAR_COMPARTMENT:
            continue

        formula = metabolite.formula or ""
        carbon_atoms = count_carbons(formula)

        if carbon_atoms == 0:
            skipped.append({
                "Reaction_ID": reaction.id,
                "Metabolite_ID": metabolite.id,
                "Metabolite_Name": metabolite.name,
                "Formula": formula,
                "Flux_mmol_gDW_h": flux,
            })
            continue

        carbon_flux = abs(flux) * carbon_atoms

        row = {
            "Reaction_ID": reaction.id,
            "Metabolite_ID": metabolite.id,
            "Metabolite_Name": metabolite.name,
            "Formula": formula,
            "Carbon_Atoms": carbon_atoms,
            "Flux_mmol_gDW_h": flux,
            "Carbon_Flux_Cmmol_gDW_h": carbon_flux,
        }

        if flux < 0:
            carbon_inputs.append(row)
        else:
            carbon_outputs.append(row)

    df_inputs = make_carbon_dataframe(carbon_inputs, "Percent_of_Total_C_Input")
    df_outputs = make_carbon_dataframe(carbon_outputs, "Percent_of_Total_C_Output")
    df_skipped = pd.DataFrame(skipped)

    return df_inputs, df_outputs, df_skipped


def make_carbon_dataframe(rows, percent_column):
    # Convert carbon-flow rows into a sorted DataFrame.
    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    df = df.sort_values("Carbon_Flux_Cmmol_gDW_h", ascending=False).reset_index(drop=True)

    total = df["Carbon_Flux_Cmmol_gDW_h"].sum()
    df[percent_column] = 100 * df["Carbon_Flux_Cmmol_gDW_h"] / total

    return df


def calculate_biomass_carbon_drain(model, solution):
    # Calculate carbon stored in biomass.
    biomass_reaction = model.reactions.get_by_id(BIOMASS_ID)
    biomass_flux = solution.fluxes[BIOMASS_ID]

    net_carbon_stoichiometry = 0

    for metabolite, coefficient in biomass_reaction.metabolites.items():
        carbon_atoms = count_carbons(metabolite.formula)
        net_carbon_stoichiometry += coefficient * carbon_atoms

    carbon_to_biomass = -net_carbon_stoichiometry * biomass_flux

    return biomass_flux, net_carbon_stoichiometry, carbon_to_biomass


def dataframe_total(df):
    # Return total carbon flux from a carbon-flow DataFrame.
    if df.empty:
        return 0

    return df["Carbon_Flux_Cmmol_gDW_h"].sum()


# Main analysis

def run_carbon_balance(model_name, model_path):
    # Load model, run pFBA, and return all carbon balance results.
    print(f"\n{model_name}: loading.....")
    model = cobra.io.read_sbml_model(str(model_path))
    model.objective = BIOMASS_ID

    print(f"{model_name}: running pFBA.....")
    solution = pfba(model)

    print(f"{model_name}: status={solution.status}  biomass={solution.fluxes[BIOMASS_ID]:.6f} h^-1")

    df_inputs, df_outputs, df_skipped = collect_exchange_carbon_fluxes(model, solution)

    total_carbon_input = dataframe_total(df_inputs)
    total_carbon_output = dataframe_total(df_outputs)

    biomass_flux, _, carbon_to_biomass = calculate_biomass_carbon_drain(model, solution)

    residual = total_carbon_input - total_carbon_output - carbon_to_biomass

    print(f"\n{model_name}: carbon balance summary")
    print(f"  Total carbon input:      {total_carbon_input:.6g} C-mmol/gDW/h")
    print(f"  Total carbon output:     {total_carbon_output:.6g} C-mmol/gDW/h")
    print(f"  Carbon to biomass:       {carbon_to_biomass:.6g} C-mmol/gDW/h")
    print(f"  Residual carbon balance: {residual:.6g} C-mmol/gDW/h")
    print(f"  C uptake exchanges:      {len(df_inputs)}")
    print(f"  C secretion exchanges:   {len(df_outputs)}")
    print(f"  Skipped exchanges:       {len(df_skipped)}")

    summary = pd.DataFrame({
        "Metric": [
            "Model",
            "Biomass reaction",
            "Biomass flux (h^-1)",
            "Total carbon input (C-mmol/gDW/h)",
            "Total carbon output (C-mmol/gDW/h)",
            "Carbon to biomass (C-mmol/gDW/h)",
            "Residual carbon balance (C-mmol/gDW/h)",
            "C-bearing uptake exchanges",
            "C-bearing secretion exchanges",
            "Skipped exchanges",
        ],
        "Value": [
            model_name,
            BIOMASS_ID,
            biomass_flux,
            total_carbon_input,
            total_carbon_output,
            carbon_to_biomass,
            residual,
            len(df_inputs),
            len(df_outputs),
            len(df_skipped),
        ],
    })

    return summary, df_inputs, df_outputs, df_skipped


OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)

with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
    for model_name, model_path in MODELS.items():
        summary, df_inputs, df_outputs, df_skipped = run_carbon_balance(model_name, model_path)

        # Each model gets its own set of sheets.
        summary.to_excel(writer, sheet_name=f"{model_name}_Summary", index=False)
        df_inputs.to_excel(writer, sheet_name=f"{model_name}_C_Inputs", index=False)
        df_outputs.to_excel(writer, sheet_name=f"{model_name}_C_Outputs", index=False)
        df_skipped.to_excel(writer, sheet_name=f"{model_name}_Exchange_Skipped", index=False)

print(f"\nWriting results to {OUTPUT_XLSX}.....")
print("Analysis done.")