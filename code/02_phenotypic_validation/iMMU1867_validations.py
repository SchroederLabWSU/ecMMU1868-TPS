'''
Phenotypic validation tests for iMMU1867.

This evaluates whether the curated model (iMMU1867) reproduces
physiologically consistent metabolic behavior.

Validation tests:
1. Biomass production using different carbon substrates
2. Respiratory quotient (RQ) during substrate oxidation
3. Biomass response under oxygen limitation
4. Single-nutrient omission analysis (from the defined minimal medium)

Output validation_results_iMMU1867.xlsx in results/02_phenotypic_validation path/:
'''

import numpy as np
import pandas as pd
from pathlib import Path
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
from cobra.exceptions import OptimizationError



# All paths are relative to the repository root.
SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_PATH = SCRIPT_DIR.parents[1] / "models" / "02_curated_iMMU1867" / "iMMU1867.xml"
OUTPUT_XLSX = SCRIPT_DIR.parents[1] / "results" / "02_phenotypic_validation" / "validation_results_iMMU1867.xlsx"

BIOMASS_ID = "BIOMASS_reaction"
ATPM_ID = "ATPM"

O2_EX = "EX_o2_e"
CO2_EX = "EX_co2_e"
HCO3_EX = "EX_hco3_e"

NH4_EX = "EX_nh4_e"
UREA_EX = "EX_urea_e"


# Defined minimal medium

DEFINED_MEDIUM = {
    # Major ions and solvent
    "EX_ca2_e": (-1000, 1000),
    "EX_cl_e": (-1000, 1000),
    "EX_h2o_e": (-1000, 1000),
    "EX_h_e": (-1000, 1000),
    "EX_k_e": (-1000, 1000),
    "EX_na1_e": (-1000, 1000),

    # Primary carbon source
    "EX_glc__D_e": (-10, 1000),

    # Inorganic nutrients and gases
    "EX_nh4_e": (-10, 1000),
    "EX_o2_e": (-10, 1000),
    "EX_pi_e": (-10, 1000),
    "EX_so4_e": (-10, 1000),

    # Essential amino acids
    "EX_his__L_e": (-0.041322, 1000),
    "EX_ile__L_e": (-0.093519, 1000),
    "EX_leu__L_e": (-0.178338, 1000),
    "EX_lys__L_e": (-0.193562, 1000),
    "EX_met__L_e": (-0.065246, 1000),
    "EX_phe__L_e": (-0.137016, 1000),
    "EX_thr__L_e": (-0.102218, 1000),
    "EX_trp__L_e": (-0.004350, 1000),
    "EX_val__L_e": (-0.115267, 1000),
}


# Carbon substrates used for phenotypic validation
# Format: substrate name, exchange ID, carbon atoms, uptake rate

SUBSTRATES = [
    ("glucose", "EX_glc__D_e", 6, 10.0), # All substrates are normalized to glucose uptake of 10 mmol/gDW/h.
    ("lactate", "EX_lac__L_e", 3, 20.0),
    ("glutamine", "EX_gln__L_e", 5, 12.0),
    ("palmitate", "EX_hdca_e", 16, 3.75),
    ("beta-hydroxybutyrate", "EX_bhb_e", 4, 15.0),
]

# Theoretical complete-oxidation reference.
THEORETICAL_RQ = pd.DataFrame([
    {"substrate": "glucose",
     "reaction": "C6H12O6 + 6 O2 -> 6 CO2 + 6 H2O",
     "RQ_theoretical": 6 / 6},
    {"substrate": "lactate",
     "reaction": "C3H6O3 + 3 O2 -> 3 CO2 + 3 H2O",
     "RQ_theoretical": 3 / 3},
    {"substrate": "glutamine",
     "reaction": "C5H10N2O3 + 4.5 O2 -> 4 CO2 + CH4N2O + 3 H2O",
     "RQ_theoretical": 4 / 4.5},
    {"substrate": "palmitate",
     "reaction": "C16H32O2 + 23 O2 -> 16 CO2 + 16 H2O",
     "RQ_theoretical": 16 / 23},
    {"substrate": "beta-hydroxybutyrate",
     "reaction": "C4H8O3 + 4.5 O2 -> 4 CO2 + 4 H2O",
     "RQ_theoretical": 4 / 4.5},
])


def close_all_exchanges(model):
    # Close uptake for all exchange reactions.
    for reaction in model.exchanges:
        reaction.lower_bound = 0
        reaction.upper_bound = 1000


def apply_medium(model, medium):
    # Apply defined minimal medium to the model.
    close_all_exchanges(model)

    for reaction_id, (lb, ub) in medium.items():
        if reaction_id in model.reactions:
            reaction = model.reactions.get_by_id(reaction_id)
            reaction.lower_bound = lb
            reaction.upper_bound = ub


def set_uptake(model, reaction_id, uptake):
    # Allow uptake for one exchange reaction.
    reaction = model.reactions.get_by_id(reaction_id)
    reaction.lower_bound = -abs(uptake)
    reaction.upper_bound = 1000


def block_exchange(model, reaction_id):
    # Block an exchange reaction.
    reaction = model.reactions.get_by_id(reaction_id)
    reaction.lower_bound = 0
    reaction.upper_bound = 0


def solve_fba(model):
    # Run FBA and return the solution.
    solution = model.optimize()

    if solution.status != "optimal":
        raise OptimizationError(solution.status)

    return solution



# Load iMMU1867 model

print("Model loading.....")

model = read_sbml_model(MODEL_PATH)
model.objective = BIOMASS_ID

print(f"Reactions:   {len(model.reactions)}")
print(f"Metabolites: {len(model.metabolites)}")
print(f"Genes:       {len(model.genes)}")


# TEST 1: Evaluate biomass production using different carbon substrates
# -----------------------------------------------------------------
print("\nTEST 1: Substrate-supported biomass")

growth_results = []

for substrate_name, exchange_id, carbons, uptake in SUBSTRATES:
    with model as m:
        apply_medium(m, DEFINED_MEDIUM)

        for _, other_exchange, _, _ in SUBSTRATES:
            if other_exchange in m.reactions:
                block_exchange(m, other_exchange)

        if exchange_id in m.reactions:
            set_uptake(m, exchange_id, uptake)

        try:
            solution = solve_fba(m)
            biomass = solution.objective_value
            status = solution.status
        except Exception:
            biomass = np.nan
            status = "failed"

        growth_results.append({
            "substrate": substrate_name,
            "exchange_id": exchange_id,
            "carbon_atoms": carbons,
            "uptake_rate": uptake,
            "biomass_h^-1": biomass,
            "status": status,
        })

df_growth = pd.DataFrame(growth_results)
print(df_growth)



# TEST 2: Evaluate substrate oxidation using respiratory quotient
#-----------------------------------------------------------------

print("\nTEST 2: Respiratory quotient")

rq_results = []

for substrate_name, exchange_id, carbons, uptake in SUBSTRATES:
    with model as m:
        close_all_exchanges(m)

        # Open inorganic exchanges needed for oxidation.
        inorganic_exchanges = [
            "EX_o2_e",
            "EX_h2o_e",
            "EX_h_e",
            "EX_na1_e",
            "EX_k_e",
            "EX_cl_e",
            "EX_ca2_e",
            "EX_nh4_e",
            "EX_pi_e",
            "EX_so4_e",
        ]

        for exchange in inorganic_exchanges:
            if exchange in m.reactions:
                m.reactions.get_by_id(exchange).lower_bound = -1000

        # Allow CO2 and bicarbonate (HCO3) secretion because some pathways can release carbon as both, 
        # and some pathways may have one but not the other. 
        for exchange in [CO2_EX, HCO3_EX]:
            if exchange in m.reactions:
                reaction = m.reactions.get_by_id(exchange)
                reaction.lower_bound = 0
                reaction.upper_bound = 1000

        # Open the selected carbon substrate.
        if exchange_id in m.reactions:
            set_uptake(m, exchange_id, uptake)

        # Glutamine-specific nitrogen handling.
        # Blocks free ammonium secretion and permits urea secretion, 
        # reflecting physiological nitrogen disposal during glutamine catabolism.
        if substrate_name == "glutamine":
            if NH4_EX in m.reactions:
                m.reactions.get_by_id(NH4_EX).upper_bound = 0

            if UREA_EX in m.reactions:
                urea = m.reactions.get_by_id(UREA_EX)
                urea.lower_bound = 0
                urea.upper_bound = 1000

        # Block biomass so substrate oxidation drives ATP maintenance.
        biomass_reaction = m.reactions.get_by_id(BIOMASS_ID)
        biomass_reaction.lower_bound = 0
        biomass_reaction.upper_bound = 0

        # Maximize ATP maintenance.
        m.objective = ATPM_ID
        if ATPM_ID in m.reactions:
            m.reactions.get_by_id(ATPM_ID).upper_bound = 1000

        try:
            solution = pfba(m)

            o2_flux = abs(solution.fluxes.get(O2_EX, np.nan))
            co2_flux = solution.fluxes.get(CO2_EX, 0)
            hco3_flux = solution.fluxes.get(HCO3_EX, 0)

            nh4_flux = solution.fluxes.get(NH4_EX, 0) if NH4_EX in m.reactions else 0
            urea_flux = solution.fluxes.get(UREA_EX, 0) if UREA_EX in m.reactions else 0

            total_co2 = co2_flux + hco3_flux
            rq = total_co2 / o2_flux if o2_flux > 0 else np.nan
            status = solution.status

        except Exception:
            o2_flux = np.nan
            co2_flux = np.nan
            hco3_flux = np.nan
            nh4_flux = np.nan
            urea_flux = np.nan
            total_co2 = np.nan
            rq = np.nan
            status = "failed"

        rq_results.append({
            "substrate": substrate_name,
            "exchange_id": exchange_id,
            "carbon_atoms": carbons,
            "uptake_rate": uptake,
            "O2_flux": o2_flux,
            "CO2_flux": co2_flux,
            "HCO3_flux": hco3_flux,
            "CO2_plus_HCO3_flux": total_co2,
            "RQ_model": rq,
            "status": status,
        })

df_rq = pd.DataFrame(rq_results)

print(df_rq)


# TEST 3: Evaluate biomass response under oxygen limitation
# ---------------------------------------------------

print("\nTEST 3: Hypoxia sweep")

O2_LIMITS = [10, 5, 2, 1, 0.5, 0.2, 0]

hypoxia_results = []

for substrate_name, exchange_id, carbons, uptake in SUBSTRATES:
    for o2_limit in O2_LIMITS:
        with model as m:
            apply_medium(m, DEFINED_MEDIUM)

            for _, other_exchange, _, _ in SUBSTRATES:
                if other_exchange in m.reactions:
                    block_exchange(m, other_exchange)

            if exchange_id in m.reactions:
                set_uptake(m, exchange_id, uptake)

            if O2_EX in m.reactions:
                oxygen = m.reactions.get_by_id(O2_EX)
                oxygen.lower_bound = -o2_limit

            try:
                solution = solve_fba(m)
                biomass = solution.objective_value
                status = solution.status
            except Exception:
                biomass = np.nan
                status = "failed"

            hypoxia_results.append({
                "substrate": substrate_name,
                "exchange_id": exchange_id,
                "O2_limit": o2_limit,
                "biomass_h^-1": biomass,
                "status": status,
            })

df_hypoxia = pd.DataFrame(hypoxia_results)
print(df_hypoxia.head())


# TEST 4: Identify essential nutrients using single-nutrient omission
# -------------------------------------------------------------------

print("\nTEST 4: Nutrient omission")

with model as m:
    apply_medium(m, DEFINED_MEDIUM)
    baseline_solution = solve_fba(m)
    baseline_biomass = baseline_solution.objective_value

omission_results = []

for reaction_id in DEFINED_MEDIUM:
    with model as m:
        apply_medium(m, DEFINED_MEDIUM)

        if reaction_id in m.reactions:
            block_exchange(m, reaction_id)

        try:
            solution = solve_fba(m)
            biomass = solution.objective_value
            fraction = biomass / baseline_biomass
            status = solution.status
        except Exception:
            biomass = np.nan
            fraction = 0
            status = "failed"

        omission_results.append({
            "omitted_exchange": reaction_id,
            "biomass_h^-1": biomass,
            "fraction_of_baseline": fraction,
            "status": status,
        })

df_omission = pd.DataFrame(omission_results)
df_omission = df_omission.sort_values("fraction_of_baseline")

print(df_omission)


# Export results to Excel

print(f"\nWriting results to {OUTPUT_XLSX}")
OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)

with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
    df_growth.to_excel(writer, sheet_name="substrate_growth", index=False)
    df_rq.to_excel(writer, sheet_name="oxidation_RQ", index=False)
    THEORETICAL_RQ.to_excel(writer, sheet_name="theoretical_RQ", index=False)
    df_hypoxia.to_excel(writer, sheet_name="hypoxia_sweep", index=False)
    df_omission.to_excel(writer, sheet_name="nutrient_omission", index=False)

print("Analysis done.")