"""
Step 2b — neuron proteomics replicate median.

Takes the log2 LFQ intensities for the three DIV15 (15 days in vitro; representing mature neurons)
neuron replicates from the Sharma et al. (2015) brain proteome and computes the per-gene median.

Input:  data/proteomics/EMS116660_proteomics_edited.xlsx
Output: results/neurons_LFQ_with_median.xlsx
"""

from pathlib import Path
import pandas as pd

SCRIPT_DIR  = Path(__file__).resolve().parent
PROTEOMICS_PATH   = SCRIPT_DIR.parents[1] / "data" / "proteomics" / "EMS116660_proteomics_edited.xlsx"
OUTPUT_XLSX = SCRIPT_DIR.parents[1] / "results" / "04_cell_specific_iMAT" / "neurons_LFQ_with_median.xlsx"

df = pd.read_excel(PROTEOMICS_PATH)

# Column names must match the proteomics sheet.
gene_col = "Gene names"
neuron_cols = ["Neurons div15 1", "Neurons div15 2", "Neurons div15 3"]

out_neuron = df[[gene_col] + neuron_cols].copy()
out_neuron["neuron_div15_log2_median"] = out_neuron[neuron_cols].median(axis=1, skipna=True)
out_neuron = out_neuron.rename(columns={gene_col: "gene_name"})

OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
out_neuron.to_excel(OUTPUT_XLSX, index=False)
print(f"Saved: {OUTPUT_XLSX}")
