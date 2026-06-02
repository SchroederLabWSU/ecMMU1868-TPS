"""
Step 2a — astrocyte proteomics replicate median.

Takes the log2 LFQ intensities for the three astrocyte replicates from the
Sharma et al. (2015) brain proteome and computes the per-gene median.

Input:  data/proteomics/EMS116660_proteomics_edited.xlsx
Output: results/04_cell_specific_iMAT/astrocytes_LFQ_with_median.xlsx
"""

from pathlib import Path
import pandas as pd

SCRIPT_DIR  = Path(__file__).resolve().parent
PROTEOMICS_PATH   = SCRIPT_DIR.parents[1] / "data" / "proteomics" / "EMS116660_proteomics_edited.xlsx"
OUTPUT_XLSX = SCRIPT_DIR.parents[1] / "results" / "04_cell_specific_iMAT" / "astrocytes_LFQ_with_median.xlsx"

df = pd.read_excel(PROTEOMICS_PATH)

# Column names must match the proteomics sheet.
gene_col = "Gene names"
astro_cols = ["Astrocytes 1", "Astrocytes 2", "Astrocytes 3"]

out_astro = df[[gene_col] + astro_cols].copy()
out_astro["astrocyte_log2_median"] = out_astro[astro_cols].median(axis=1, skipna=True)
out_astro = out_astro.rename(columns={gene_col: "gene_name"})

OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
out_astro.to_excel(OUTPUT_XLSX, index=False)
print(f"Saved: {OUTPUT_XLSX}")
