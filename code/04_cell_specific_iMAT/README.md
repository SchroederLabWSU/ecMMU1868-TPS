# Cell-specific reconstruction using integrative Metabolic Analysis Tool (iMAT)

Derive astrocyte-specific (**iMMU1867a**) and neuron-specific (**iMMU1867n**)
models from iMMU1867 by integrating cell-type proteomics with iMAT. Also computes
MEMOTE quality scores and active-reaction overlap.

Sub-steps (mirroring the proteomics → iMAT workflow):
1. Map model gene names → NCBI gene IDs (`../../data/annotations/`).
2. Compute replicate-median protein abundances (`../../data/proteomics/`).
3. Merge model genes with proteomics medians; report % mapped.
4. Run the iMAT workflow to extract cell-specific models.

- **Input:**  `../../models/01_curated_iMMU1867/iMMU1867.xml`,
  `../../data/proteomics/`, `../../data/annotations/`
- **Output:** `../../models/02_cell_specific_iMAT/iMMU1867a.xml`,
  `iMMU1867n.xml`; panels for **Fig 4**
- **Entry point:** _add your iMAT workflow here_ (e.g. `run_imat.py` / notebook)
