"""
Active Reaction Overlap — Venn Diagram + Excel Export 

Compares the active reactions of three models:
  - iMMU1867   : curated generic Mus musculus model
  - iMMU1867a  : iMAT astrocyte-specific reconstruction
  - iMMU1867n  : iMAT neuron-specific reconstruction

Runs pFBA on each, plots a Venn diagram of their active-reaction overlap, and
saves all active reactions to Excel.

Outputs (results/04_cell_specific_iMAT/):
  - active_reaction_overlap.png
  - active_reaction_overlap.xlsx
"""

from pathlib import Path
import cobra
from cobra.flux_analysis import pfba
import matplotlib.pyplot as plt
from matplotlib_venn import venn3_unweighted
import pandas as pd


# All paths are relative to the repository root.
SCRIPT_DIR  = Path(__file__).resolve().parent
MODELS_DIR  = SCRIPT_DIR.parents[1] / "models"
RESULTS_DIR = SCRIPT_DIR.parents[1] / "results" / "04_cell_specific_iMAT"

MODEL_PATHS = {
    "iMMU1867": MODELS_DIR / "02_curated_iMMU1867"  / "iMMU1867.xml",
    "iMMU1867a": MODELS_DIR / "03_cell_specific_iMAT" / "iMMU1867a.xml",
    "iMMU1867n": MODELS_DIR / "03_cell_specific_iMAT" / "iMMU1867n.xml",
}
BIOMASS_RXN_ID = "BIOMASS_reaction"
FLUX_TOL = 1e-6 


def get_active_reactions(model_path, label):
    # Load model, run pFBA and return (active_id_set, rxn_fluxes, rxn_names).
    model = cobra.io.read_sbml_model(str(model_path))
    model.objective = BIOMASS_RXN_ID

    sol = pfba(model)
    active = set(sol.fluxes.index[sol.fluxes.abs() > FLUX_TOL])
    names = {rxn.id: rxn.name for rxn in model.reactions}
    print(f"  {label} - Biomass : {sol.fluxes[BIOMASS_RXN_ID]:.6f} 1/h | active reactions: {len(active)}")
    return active, sol.fluxes, names


# Run pFBA on each model and get the active reactions, fluxes, and names.
active, fluxes, names = {}, {}, {}
for label, path in MODEL_PATHS.items():
    active[label], fluxes[label], names[label] = get_active_reactions(path, label)

# All three models share the same reaction namespace.
rxn_names = {}
for name_map in names.values():
    rxn_names.update(name_map)

A = active["iMMU1867"]
B = active["iMMU1867a"]
C = active["iMMU1867n"]

# The Venn regions
regions = {
    "iMMU1867_only": A - (B | C),
    "iMMU1867a_only": B - (A | C),
    "iMMU1867n_only": C - (A | B),
    "iMMU1867_iMMU1867a": (A & B) - C,
    "iMMU1867_iMMU1867n": (A & C) - B,
    "iMMU1867a_iMMU1867n":  (B & C) - A,
    "All_three": A & B & C,
}
overlap_counts = {k: len(v) for k, v in regions.items()}

print("\n--- Active Reaction Overlap ---")
for k, v in overlap_counts.items():
    print(f"  {k:<22s}: {v}")

# Write to Excel
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
excel_path = RESULTS_DIR / "active_reaction_overlap.xlsx"

with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
    # Summary counts
    pd.DataFrame(list(overlap_counts.items()),
                 columns=["Overlap_Group", "Reaction_Count"]
                 ).to_excel(writer, sheet_name="Summary", index=False)
    # One sheet per Venn region
    for region, rxn_set in regions.items():
        ids = sorted(rxn_set)
        pd.DataFrame({
            "Reaction_ID":   ids,
            "Reaction_Name": [rxn_names.get(r, "") for r in ids],
        }).to_excel(writer, sheet_name=region[:31], index=False)

print(f"\nExcel workbook saved: {excel_path}")

# Venn diagram
# venn3 subset order: (A_only, B_only, AB, C_only, AC, BC, ABC)
subsets = (
    overlap_counts["iMMU1867_only"],
    overlap_counts["iMMU1867a_only"],
    overlap_counts["iMMU1867_iMMU1867a"],
    overlap_counts["iMMU1867n_only"],
    overlap_counts["iMMU1867_iMMU1867n"],
    overlap_counts["iMMU1867a_iMMU1867n"],
    overlap_counts["All_three"],
)

# Plot and save the Venn diagram.
fig, ax = plt.subplots(figsize=(9, 6))
v = venn3_unweighted(subsets=subsets,
     set_labels=("iMMU1867", "iMMU1867a", "iMMU1867n"), ax=ax)
for txt in v.set_labels:
    if txt:
        txt.set_fontsize(14)
for txt in v.subset_labels:
    if txt:
        txt.set_fontsize(12)
ax.set_title("Active Reaction Overlap (pFBA)", fontsize=14, pad=16)

fig_path = RESULTS_DIR / "active_reaction_overlap.png"
plt.tight_layout()
plt.savefig(fig_path, dpi=150, bbox_inches="tight")
print(f"Venn diagram saved:  {fig_path}")
