"""
Model statistics across reconstruction stages.

For each model reports reaction, metabolite, gene, and compartment counts
and the optimal biomass production. The model consists of the curated generic model,
the iMAT cell-specific models, and the tripartite-expanded models.

Output: results/05_tripartite_synapse_expansion/model_statistics.xlsx
"""

from pathlib import Path
import cobra
import pandas as pd

# All paths are relative to the repository root.
SCRIPT_DIR  = Path(__file__).resolve().parent
MODELS_DIR  = SCRIPT_DIR.parents[1] / "models"
OUTPUT_XLSX = SCRIPT_DIR.parents[1] / "results" / "05_tripartite_synapse_expansion" / "model_statistics.xlsx"

MODELS = {
    "iMMU1867": MODELS_DIR / "02_curated_iMMU1867"   / "iMMU1867.xml",  # curated generic model
    "iMMU1867a": MODELS_DIR / "03_cell_specific_iMAT"  / "iMMU1867a.xml",  # iMAT astrocyte
    "iMMU1867n": MODELS_DIR / "03_cell_specific_iMAT"  / "iMMU1867n.xml",  # iMAT neuron
    "iMMU1868a-pap": MODELS_DIR / "04_tripartite_expanded" / "iMMU1868a-pap.xml",  # astrocyte + PAP
    "iMMU1868n-pre": MODELS_DIR / "04_tripartite_expanded" / "iMMU1868n-pre.xml",  # pre-synaptic neuron
    "iMMU1868n-post": MODELS_DIR / "04_tripartite_expanded" / "iMMU1868n-post.xml", # post-synaptic neuron
}

BIOMASS_RXN_ID = "BIOMASS_reaction"

rows = []
for name, path in MODELS.items():
    print(f"Loading {name} .....")
    model = cobra.io.read_sbml_model(str(path))
    model.objective = BIOMASS_RXN_ID

    sol = model.optimize()
    biomass = sol.objective_value

    rows.append({
        "Model": name,
        "Reactions": len(model.reactions),
        "Metabolites": len(model.metabolites),
        "Genes": len(model.genes),
        "Compartments":  len(model.compartments),
        "Biomass (1/h)":  round(biomass, 4),
    })

df = pd.DataFrame(rows)
print("\n", df.to_string(index=False))

OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
df.to_excel(OUTPUT_XLSX, index=False)
print(f"\nSaved: {OUTPUT_XLSX}")
