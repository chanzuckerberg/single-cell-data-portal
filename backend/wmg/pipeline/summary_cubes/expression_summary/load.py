import logging

import numpy as np
import pandas as pd

from backend.wmg.pipeline.summary_cubes.cell_count import remove_accents, return_dataset_dict_w_publications

logger: logging.Logger = logging.getLogger(__name__)


def build_in_mem_cube(
    *,
    gene_ids: pd.DataFrame,
    cube_index: pd.DataFrame,
    other_cube_attrs: list,
    cube_sum: np.ndarray,
    cube_nnz: np.ndarray,
):
    """
    Build the cube in memory, calculating the gene expression value for each combination of attributes
    """
    logger.info("Building in-mem cube")

    # Count total values so we can allocate buffers once
    total_vals = 0
    for cube_idx in cube_index.cube_idx.values:
        mask = cube_nnz[cube_idx] != 0
        total_vals += np.count_nonzero(mask)

    # allocate buffers
    dims = [
        np.empty((total_vals,), dtype=object),
        np.empty((total_vals,), dtype=object),
        np.empty((total_vals,), dtype=object),
    ]
    vals = {
        "sum": np.empty((total_vals,)),
        "nnz": np.empty((total_vals,), dtype=np.uint64),
        **{k: np.empty((total_vals,), dtype=object) for k in other_cube_attrs},
    }

    # populate buffers
    idx = 0
    dataset_dict = return_dataset_dict_w_publications()
    for grp in cube_index.to_records():
        (
            tissue_ontology_term_id,
            organism_ontology_term_id,
            *attr_values,
            _,
            cube_idx,
        ) = grp.tolist()
        mask = cube_nnz[cube_idx] != 0
        n_vals = np.count_nonzero(mask)
        if n_vals == 0:  # Used to maintain sparsity
            continue

        logger.debug(grp)

        dims[0][idx : idx + n_vals] = gene_ids.gene_ontology_term_id.values[mask]
        dims[1][idx : idx + n_vals] = tissue_ontology_term_id
        dims[2][idx : idx + n_vals] = organism_ontology_term_id

        vals["sum"][idx : idx + n_vals] = cube_sum[cube_idx, mask]
        vals["nnz"][idx : idx + n_vals] = cube_nnz[cube_idx, mask]

        dataset_index = 0
        for i, k in enumerate(other_cube_attrs):
            if k == "dataset_id":
                dataset_index = i
            if k != "publication_citation":
                vals[k][idx : idx + n_vals] = attr_values[i]

        vals["publication_citation"][idx : idx + n_vals] = remove_accents(
            dataset_dict.get(attr_values[dataset_index].decode(), "No Publication")
        )

        idx += n_vals

    return dims, vals


def build_in_mem_cube_default(
    *,
    gene_ids: pd.DataFrame,
    cube_index: pd.DataFrame,
    other_cube_attrs: list,
    cube_sum: np.ndarray,
    cube_nnz: np.ndarray,
    cube_sqsum: np.ndarray,
):
    """
    Build the cube in memory, calculating the gene expression value for each combination of attributes
    """
    logger.info("Building in-mem default cube")

    # Count total values so we can allocate buffers once
    total_vals = 0
    for cube_idx in cube_index.cube_idx.values:
        mask = cube_nnz[cube_idx] != 0
        total_vals += np.count_nonzero(mask)

    # allocate buffers
    dims = [
        np.empty((total_vals,), dtype=object),
        np.empty((total_vals,), dtype=object),
        np.empty((total_vals,), dtype=object),
    ]
    vals = {
        "sum": np.empty((total_vals,)),
        "nnz": np.empty((total_vals,), dtype=np.uint64),
        "sqsum": np.empty((total_vals,)),
        **{k: np.empty((total_vals,), dtype=object) for k in other_cube_attrs},
    }

    # populate buffers
    idx = 0

    for grp in cube_index.to_records():
        (
            tissue_ontology_term_id,
            organism_ontology_term_id,
            *attr_values,
            _,
            cube_idx,
        ) = grp.tolist()
        mask = cube_nnz[cube_idx] != 0
        n_vals = np.count_nonzero(mask)
        if n_vals == 0:  # Used to maintain sparsity
            continue

        logger.debug(grp)

        dims[0][idx : idx + n_vals] = gene_ids.gene_ontology_term_id.values[mask]
        dims[1][idx : idx + n_vals] = tissue_ontology_term_id
        dims[2][idx : idx + n_vals] = organism_ontology_term_id

        vals["sum"][idx : idx + n_vals] = cube_sum[cube_idx, mask]
        vals["nnz"][idx : idx + n_vals] = cube_nnz[cube_idx, mask]
        vals["sqsum"][idx : idx + n_vals] = cube_sqsum[cube_idx, mask]

        for i, k in enumerate(other_cube_attrs):
            vals[k][idx : idx + n_vals] = attr_values[i]

        idx += n_vals

    return dims, vals
