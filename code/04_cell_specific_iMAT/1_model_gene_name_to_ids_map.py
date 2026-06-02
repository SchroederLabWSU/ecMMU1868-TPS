"""
Step 1 — map model gene symbols to NCBI Gene IDs.

Parses the curated model (iMMU1867) gene symbols directly to NCBI-Gene-ID 
table for proteomics integration.

Input:  models/02_curated_iMMU1867/iMMU1867.xml
Output: results/04_cell_specific_iMAT/model_gene_name_to_ncbi_id.xlsx
"""

import re
import pandas as pd
from pathlib import Path


# All paths are relative to the repository root.
SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_PATH = SCRIPT_DIR.parents[1] / "models" / "02_curated_iMMU1867" / "iMMU1867.xml"
OUTPUT_XLSX = SCRIPT_DIR.parents[1] / "results" / "04_cell_specific_iMAT" / "model_gene_name_to_ncbi_id.xlsx"

# Read the raw model SBML so we can parse the gene annotations.
text = MODEL_PATH.read_text(encoding="utf-8", errors="ignore")

# Extract each complete fbc:geneProduct block (gene id, name, RDF annotation).
blocks = re.findall(
    r"<fbc:geneProduct\b.*?</fbc:geneProduct>",
    text,
    flags=re.IGNORECASE | re.DOTALL,
)

rows = []
for b in blocks:
    m_name = re.search(r'\bfbc:name="([^"]+)"', b)  # gene symbol
    m_ncbi = re.search(r"identifiers\.org/ncbigene/(\d+)", b)  # NCBI Gene ID
    if m_name and m_ncbi:
        rows.append((m_name.group(1), m_ncbi.group(1)))

df = pd.DataFrame(rows, columns=["gene_name", "gene_id"])
df = df.drop_duplicates(subset=["gene_id"]).sort_values(["gene_name", "gene_id"])

OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
df.to_excel(OUTPUT_XLSX, index=False)
print(f"Saved mapping with {len(df)} genes {OUTPUT_XLSX}")
