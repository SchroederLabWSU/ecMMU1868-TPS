"""
Loads a MATLAB v7.3 (HDF5) COBRA-style community model into a COBRApy Model.

We needed this because:
  - scipy.io.loadmat does not support MATLAB v7.3 (.mat) files.
  - The tripartite synapse community models could not be exported to SBML via writeCbModel.
"""

from __future__ import annotations

import h5py
import numpy as np
import cobra
from scipy.sparse import csc_matrix


def _decode_str(arr) -> str:
    # Decode a MATLAB string stored as a uint16 HDF5 array.
    if hasattr(arr, "dtype") and arr.dtype.kind in ("u", "i"):
        chars = np.asarray(arr).flatten()
        return "".join(chr(int(c)) for c in chars if int(c) != 0)
    if isinstance(arr, bytes):
        return arr.decode("utf-8", errors="replace")
    return str(arr)


def _read_cell_of_str(hdf5_file: h5py.File, dataset) -> list[str]:
    # Read a MATLAB cell array of strings from HDF5 object references.
    result: list[str] = []
    if isinstance(dataset, h5py.Dataset) and dataset.dtype == h5py.ref_dtype:
        for ref in np.asarray(dataset).flatten():
            if not ref:
                result.append("")
            else:
                result.append(_decode_str(hdf5_file[ref][()]))
    else:
        result = [_decode_str(dataset[()])]
    return result


def _read_sparse_S(group) -> csc_matrix:
    # Read a MATLAB sparse matrix stored in HDF5 CSC format.
    # MATLAB stores S as (mets x rxns) but the HDF5 layout is transposed.
    nrows = group.attrs["MATLAB_sparse"]
    data  = np.asarray(group["data"]).flatten()
    ir    = np.asarray(group["ir"]).flatten()
    jc    = np.asarray(group["jc"]).flatten()
    ncols = jc.size - 1
    return csc_matrix((data, ir, jc), shape=(int(nrows), int(ncols)))


def load_mat73_cobra(mat_path: str, var: str = "modelCom") -> cobra.Model:
    # Build a COBRApy Model from a MATLAB v7.3 MAT file containing a COBRA struct.
    hdf5_file = h5py.File(mat_path, "r")
    group     = hdf5_file[var]

    rxn_ids  = _read_cell_of_str(hdf5_file, group["rxns"])
    met_ids  = _read_cell_of_str(hdf5_file, group["mets"])
    lb       = np.asarray(group["lb"]).flatten()
    ub       = np.asarray(group["ub"]).flatten()
    obj_coef = np.asarray(group["c"]).flatten()
    b        = np.asarray(group["b"]).flatten() if "b" in group else np.zeros(len(met_ids))

    rxn_names = _read_cell_of_str(hdf5_file, group["rxnNames"]) if "rxnNames" in group else rxn_ids
    met_names = _read_cell_of_str(hdf5_file, group["metNames"]) if "metNames" in group else met_ids
    gr_rules  = _read_cell_of_str(hdf5_file, group["grRules"])  if "grRules"  in group else [""] * len(rxn_ids)
    subsystems= _read_cell_of_str(hdf5_file, group["subSystems"])if "subSystems" in group else [""] * len(rxn_ids)

    S = _read_sparse_S(group["S"])
    if S.shape[0] != len(met_ids) or S.shape[1] != len(rxn_ids):
        S = S.T.tocsc()

    model = cobra.Model("community_ec")

    # Add metabolites
    metabolites = []
    for i, met_id in enumerate(met_ids):
        met = cobra.Metabolite(id=met_id,
                               name=met_names[i] if i < len(met_names) else met_id)
        metabolites.append(met)
    model.add_metabolites(metabolites)

    # Add reactions
    reactions = []
    for j, rxn_id in enumerate(rxn_ids):
        rxn = cobra.Reaction(id=rxn_id,
                             name=rxn_names[j] if j < len(rxn_names) else rxn_id)
        rxn.lower_bound = float(lb[j])
        rxn.upper_bound = float(ub[j])
        if j < len(subsystems) and subsystems[j]:
            rxn.subsystem = subsystems[j]
        if j < len(gr_rules) and gr_rules[j]:
            rxn.gene_reaction_rule = gr_rules[j]
        reactions.append(rxn)
    model.add_reactions(reactions)

    # Apply stoichiometry from S (mets x rxns)
    S_csc = S.tocsc()
    for j in range(S_csc.shape[1]):
        start, end = S_csc.indptr[j], S_csc.indptr[j + 1]
        if start == end:
            continue
        coefficients = {metabolites[int(S_csc.indices[k])]: float(S_csc.data[k])
                        for k in range(start, end)}
        reactions[j].add_metabolites(coefficients)

    # Set objective from c vector
    obj_rxns = np.where(np.abs(obj_coef) > 0)[0]
    if obj_rxns.size:
        model.objective = {reactions[int(i)]: float(obj_coef[int(i)]) for i in obj_rxns}

    hdf5_file.close()
    return model


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python mat73_to_cobra.py <path_to_mat_file>")
        sys.exit(1)
    mat_path = sys.argv[1]
    model = load_mat73_cobra(mat_path)
    print(f"Loaded: {len(model.reactions)} rxns, {len(model.metabolites)} mets, {len(model.genes)} genes")
