"""
Model statistics and pFBA-active reactions across stoichiometric models.

Reports total reactions, metabolites, genes, compartments, pFBA-active
reactions, and maximum growth rate for each model across reconstruction
stages: curated parent, iMAT cell-specific, and tripartite expansion.

Active reactions are defined as reactions with |flux| > 1e-6 in the pFBA solution.

Output: results/05_tripartite_synapse_expansion/model_statistics_active_pFBA.xlsx
"""

from pathlib import Path
import cobra
from cobra.flux_analysis import pfba
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
MODELS_DIR = REPO_ROOT / "models"
OUTPUT_XLSX = REPO_ROOT / "results" / "05_tripartite_synapse_expansion" / "model_statistics_active_pFBA.xlsx"

MODELS = {
    "iMMU1867": MODELS_DIR / "02_curated_iMMU1867" / "iMMU1867.xml", # curated generic model
    "iMMU1867a": MODELS_DIR / "03_cell_specific_iMAT" / "iMMU1867a.xml", # iMAT astrocyte
    "iMMU1867n": MODELS_DIR / "03_cell_specific_iMAT" / "iMMU1867n.xml", # iMAT neuron
    "iMMU1868a-pap": MODELS_DIR / "04_tripartite_expanded" / "iMMU1868a-pap.xml", # astrocyte + PAP
    "iMMU1868n-pre": MODELS_DIR / "04_tripartite_expanded" / "iMMU1868n-pre.xml", # pre-synaptic neuron
    "iMMU1868n-post": MODELS_DIR / "04_tripartite_expanded" / "iMMU1868n-post.xml", # post-synaptic neuron
}

BIOMASS_ID = "BIOMASS_reaction"
ACTIVE_FLUX_THRESHOLD = 1e-6

rows = []
for name, path in MODELS.items():
    print(f"Loading {name} .....", flush=True)
    model = cobra.io.read_sbml_model(str(path))
    model.objective = BIOMASS_ID

    growth = model.optimize().objective_value
    sol = pfba(model)
    active_rxns = int((sol.fluxes.abs() > ACTIVE_FLUX_THRESHOLD).sum())

    rows.append({
        "Model": name,
        "Reactions": len(model.reactions),
        "Metabolites": len(model.metabolites),
        "Genes": len(model.genes),
        "Compartments": len(model.compartments),
        "Active reactions (pFBA)": active_rxns,
        "Maximum growth rate (h-1)": round(growth, 4),
    })
    print(f"  active rxns={active_rxns}  growth={round(growth, 4)}", flush=True)

df = pd.DataFrame(rows)
print("\n", df.to_string(index=False))

# Save to Excel.
OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
df.to_excel(OUTPUT_XLSX, index=False)
print(f"\nSaved: {OUTPUT_XLSX}")
