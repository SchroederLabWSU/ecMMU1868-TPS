"""
step 3 - merge model genes with astrocyte and neuron proteomics medians.

Creates one excel table containing each model gene, its NCBI Gene ID, and the
astrocyte and DIV15 neuron log2 LFQ median values used for iMAT
context-specific model extraction in the workflow code (iMAT_workflow.ipynb).

Inputs (from step 1, step 2a, and 2b, in results/04_cell_specific_iMAT):
  - model_gene_name_to_ncbi_id.xlsx
  - astrocytes_LFQ_with_median.xlsx
  - neurons_LFQ_with_median.xlsx
Output:
  - results/04_cell_specific_iMAT/model_gene_to_proteomics_medians.xlsx
"""

from pathlib import Path
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR.parents[1] / "results" / "04_cell_specific_iMAT"
OUTPUT_XLSX = RESULTS_DIR / "model_gene_to_proteomics_medians.xlsx"

gene_mapping_df    = pd.read_excel(RESULTS_DIR / "model_gene_name_to_ncbi_id.xlsx")
astro_df  = pd.read_excel(RESULTS_DIR / "astrocytes_LFQ_with_median.xlsx")
neuron_df = pd.read_excel(RESULTS_DIR / "neurons_LFQ_with_median.xlsx")

if "gene_name" not in gene_mapping_df.columns and "Gene names" in gene_mapping_df.columns:
    gene_mapping_df = gene_mapping_df.rename(columns={"Gene names": "gene_name"})

if "astrocyte_log2_median" not in astro_df.columns:
    for c in astro_df.columns:
        if "astro" in c.lower() and "median" in c.lower():
            astro_df = astro_df.rename(columns={c: "astrocyte_log2_median"})
            break

if "neuron_div15_log2_median" not in neuron_df.columns:
    for c in neuron_df.columns:
        if "neuron" in c.lower() and "median" in c.lower():
            neuron_df = neuron_df.rename(columns={c: "neuron_div15_log2_median"})
            break

# Drop duplicate genes.
gene_mapping_df = gene_mapping_df[["gene_name", "gene_id"]].drop_duplicates(subset=["gene_name"])
astro_df = astro_df[["gene_name", "astrocyte_log2_median"]].drop_duplicates(subset=["gene_name"])
neuron_df = neuron_df[["gene_name", "neuron_div15_log2_median"]].drop_duplicates(subset=["gene_name"])

# Keep all model genes and add proteomics medians when available.
out = (gene_mapping_df
       .merge(astro_df, on="gene_name", how="left")
       .merge(neuron_df, on="gene_name", how="left"))
out = out[["gene_name", "gene_id", "astrocyte_log2_median", "neuron_div15_log2_median"]]

OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
out.to_excel(OUTPUT_XLSX, index=False)
print(f"Saved: {OUTPUT_XLSX} with {len(out)} rows")
