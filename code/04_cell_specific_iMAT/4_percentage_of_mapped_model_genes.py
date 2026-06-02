"""
Step 4 — fraction of model genes covered by the proteomics data.

Checks how many model genes received an astrocyte or neuron median abundance
and exports a summary and the genes that differ between cell types.

Input:  results/04_cell_specific_iMAT/model_gene_to_proteomics_medians.xlsx (from step 3)
Output: results/04_cell_specific_iMAT/gene_coverage_summary.xlsx
"""

from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
MERGED_PATH = SCRIPT_DIR.parents[1] / "results" / "04_cell_specific_iMAT" / "model_gene_to_proteomics_medians.xlsx"
OUTPUT_XLSX = SCRIPT_DIR.parents[1] / "results" / "04_cell_specific_iMAT" / "gene_coverage_summary.xlsx"

df = pd.read_excel(MERGED_PATH)
total_genes = len(df)

astro_mapped  = df["astrocyte_log2_median"].notna()
neuron_mapped = df["neuron_div15_log2_median"].notna()

n_astro  = astro_mapped.sum()
n_neuron = neuron_mapped.sum()

print(f"Total model genes: {total_genes}")
print(f"Astrocyte-mapped genes: {n_astro} ({100 * n_astro / total_genes:.2f}%)")
print(f"Neuron (div15)-mapped genes: {n_neuron} ({100 * n_neuron / total_genes:.2f}%)")

# Gene groups by detection status
mapped_both = df[astro_mapped & neuron_mapped][["gene_name", "gene_id", "astrocyte_log2_median", "neuron_div15_log2_median"]]
astro_only = df[astro_mapped & ~neuron_mapped][["gene_name", "gene_id", "astrocyte_log2_median"]]
neuron_only = df[~astro_mapped & neuron_mapped][["gene_name", "gene_id", "neuron_div15_log2_median"]]
unmapped = df[~astro_mapped & ~neuron_mapped][["gene_name", "gene_id"]]

print(f"Mapped in both astrocyte and neuron: {len(mapped_both)} ({100 * len(mapped_both) / total_genes:.2f}%)")
print(f"Astrocyte-only genes: {len(astro_only)} ({100 * len(astro_only) / total_genes:.2f}%)")
print(f"Neuron-only genes: {len(neuron_only)} ({100 * len(neuron_only) / total_genes:.2f}%)")
print(f"Unmapped: {len(unmapped)} ({100 * len(unmapped) / total_genes:.2f}%)")

# Summary table
summary = pd.DataFrame({
    "Category": [
        "Total model genes",
        "Mapped in both astrocyte and neuron",
        "Astrocyte-mapped (total)",
        "Neuron-mapped (total)",
        "Astrocyte only (not in neuron)",
        "Neuron only (not in astrocyte)",
        "Unmapped in both",
    ],
    "Count": [total_genes, len(mapped_both), n_astro, n_neuron,
              len(astro_only), len(neuron_only), len(unmapped)],
    "Percentage": [
        "100%",
        f"{100 * len(mapped_both) / total_genes:.2f}%",
        f"{100 * n_astro / total_genes:.2f}%",
        f"{100 * n_neuron / total_genes:.2f}%",
        f"{100 * len(astro_only) / total_genes:.2f}%",
        f"{100 * len(neuron_only) / total_genes:.2f}%",
        f"{100 * len(unmapped) / total_genes:.2f}%",
    ],
})

# Excel export

OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
    summary.to_excel(writer, sheet_name="Summary", index=False)
    mapped_both.to_excel(writer, sheet_name="Mapped_both", index=False)
    astro_only.to_excel(writer, sheet_name="Astrocyte_only", index=False)
    neuron_only.to_excel(writer, sheet_name="Neuron_only", index=False)
    unmapped.to_excel(writer, sheet_name="Unmapped", index=False)

print(f"Saved: {OUTPUT_XLSX}")
