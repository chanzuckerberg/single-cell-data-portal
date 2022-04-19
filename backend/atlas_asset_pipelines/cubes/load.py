import logging
import time

import numpy as np
import tiledb

from backend.atlas_asset_pipelines.cubes.transform import reduce_X
from backend.wmg.data.schemas.cube_schema import cube_non_indexed_dims
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.wmg_cube import make_cube_index


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


cube_indexed_dims_no_gene_ontology = [
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
]


def load_data_into_cube(tdb_group, uri: str):
    """
    Load data from the concat_corpus into the queryable expression summary cube
    """
    ctx = create_ctx()
    start_time = time.time()
    logger.debug(f"Start loading big cube at : {uri}")

    with tiledb.open(f"{tdb_group}/var", ctx=ctx) as var:
        gene_ontology_term_ids = var.query(dims=["gene_ontology_term_id"], attrs=["var_idx"], use_arrow=False).df[:]
        gene_ontology_term_ids.sort_values(by="var_idx", inplace=True)
    n_genes = len(gene_ontology_term_ids)

    ##
    # Reduce X
    ##
    big_cube_atts = cube_indexed_dims_no_gene_ontology + cube_non_indexed_dims
    cell_labels, cube_index = make_cube_index(tdb_group, big_cube_atts)
    n_groups = len(cube_index)

    cube_sum = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_nnz = np.zeros((n_groups, n_genes), dtype=np.uint64)
    cube_min = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_max = np.zeros((n_groups, n_genes), dtype=np.float32)

    # pass 1 - sum, nnz, min, max
    reduce_X(tdb_group, start_time, cell_labels.cube_idx.values, cube_sum, cube_nnz, cube_min, cube_max)

    return build_in_mem_cube(gene_ontology_term_ids, cube_index, cube_non_indexed_dims, cube_sum, cube_nnz)


def build_in_mem_cube(gene_ids, cube_index, other_attrs, cube_sum, cube_nnz):
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
        "n_cells": np.empty((total_vals,), dtype=np.uint32),
        **{k: np.empty((total_vals,), dtype=object) for k in other_attrs},
    }

    # populate buffers
    idx = 0

    for grp in cube_index.to_records():
        (
            tissue_ontology_term_id,
            organism_ontology_term_id,
            *attr_values,
            n,
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
        vals["n_cells"][idx : idx + n_vals] = n  # wasteful

        for i, k in enumerate(other_attrs):
            vals[k][idx : idx + n_vals] = attr_values[i]

        idx += n_vals

    return dims, vals