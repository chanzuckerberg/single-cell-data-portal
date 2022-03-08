"""
TODO:
1. parallelize cube building?  See https://more-itertools.readthedocs.io/en/stable/api.html#more_itertools.chunked
2. calculate bins & mean more efficiently in cube builder
3. large cube tests

"""
import logging
import time
import concurrent.futures

import tiledb
import numpy as np
import pandas as pd

from backend.wmg.data.config import create_fast_ctx
from .compute import coo_cube_pass1_into
from ..data.cube_schema import schema, cube_indexed_dims, cube_non_indexed_dims, cube_indexed_dims_no_gene_ontology

logger = logging.getLogger(__name__)


# big_cube_attrs = [
#     "dataset_id",
#     "assay_ontology_term_id",
#     "development_stage_ontology_term_id",
#     "disease_ontology_term_id",
#     "ethnicity_ontology_term_id",
#     "sex_ontology_term_id",
# ]
# big_cube_dims = [
#     "cell_type_ontology_term_id",
#     "tissue_ontology_term_id",
#     "organism_ontology_term_id",
#     *big_cube_attrs,
# ]


# def create_empty_cube(uri: str, other_attrs=[]):
#     if tiledb.array_exists(uri):
#         tiledb.remove(uri)
#
#     filters = [tiledb.ZstdFilter(level=-22)]
#     domain = tiledb.Domain(
#         [
#             tiledb.Dim(name="feature_id", domain=None, tile=None, dtype="ascii", filters=filters),
#             tiledb.Dim(name="cell_type_ontology_term_id", domain=None, tile=None, dtype="ascii", filters=filters),
#             tiledb.Dim(name="tissue_ontology_term_id", domain=None, tile=None, dtype="ascii", filters=filters),
#             tiledb.Dim(name="organism_ontology_term_id", domain=None, tile=None, dtype="ascii", filters=filters),
#         ]
#     )
#     attrs = [
#         tiledb.Attr(name="n_cells", dtype=np.uint32, filters=filters),
#         tiledb.Attr(name="nnz", dtype=np.uint64, filters=filters),
#         tiledb.Attr(name="sum", dtype=np.float32, filters=filters),
#         *[tiledb.Attr(name=a, dtype="ascii", var=True, filters=filters) for a in other_attrs],
#     ]
#
#     tiledb.Array.create(
#         uri,
#         tiledb.ArraySchema(
#             domain=domain,
#             sparse=True,
#             allows_duplicates=True,
#             attrs=attrs,
#             cell_order="row-major",
#             tile_order="row-major",
#             capacity=10000,
#         ),
#     )


def make_cube_index(tdb_group, cube_dims):
    """
    Create index for queryable dimensions
    """
    with tiledb.open(f"{tdb_group}/obs") as obs:
        cell_labels = obs.query(use_arrow=False).df[:]
    cell_labels.sort_values(by=['obs_idx'], inplace=True, ignore_index=True)

    cell_labels = pd.DataFrame(
        # Todo previously cube dimensions did not include the feature_id (gene_ontology_id) now its  included ... breaking change?
        data={k: cell_labels[k].astype("category") for k in cube_dims},
        index=cell_labels.obs_idx,
    )

    cube_index = pd.DataFrame(cell_labels.value_counts(), columns=["n"])
    cube_index["cube_idx"] = range(len(cube_index))

    cell_labels = cell_labels.join(cube_index.cube_idx, on=cube_dims)

    # we failed to correctly create the corpus if these are false
    assert len(cell_labels.index) == cell_labels.index[-1] + 1
    assert cell_labels.index[0] == 0

    return cell_labels, cube_index


def reduce_X(tdb_group, start_time, cube_indices, *accum):
    """
    # TODO
    """
    with concurrent.futures.ThreadPoolExecutor() as tp:
        cfg = {
            "py.init_buffer_bytes": 512 * 1024 ** 2,
            "py.exact_init_buffer_bytes": "true",
        }
        with tiledb.open(f"{tdb_group}/raw", ctx=create_fast_ctx(cfg)) as X:
            iterable = X.query(return_incomplete=True, order="U", attrs=["data"])
            future = None
            for i, result in enumerate(iterable.df[:]):
                logger.info(f"reduce raw X, iter {i}, {time.time()-start_time}")
                if future is not None:
                    future.result()  # forces a wait
                future = tp.submit(
                    coo_cube_pass1_into,
                    result["data"].values,
                    result["obs_idx"].values,
                    result["var_idx"].values,
                    cube_indices,
                    *accum,
                )

        return accum


def build_in_mem_cube(gene_ids, cube_index, other_attrs, cube_sum, cube_nnz):
    """
    Build the cube in memory, calculating the
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
        if n_vals == 0: # Todo what is the purpose of this?
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


def load_data_into_cube(tdb_group, uri: str):
    """
    Load data from the concatenated corpus into the queryable cube
    """
    ctx = create_fast_ctx()
    start_time = time.time()
    logger.debug(f"Start loading big cube at : {uri}")

    with tiledb.open(f"{tdb_group}/var", ctx=ctx) as var:
        gene_ontology_term_ids = var.query(dims=["gene_ontology_term_id"], attrs=[], use_arrow=False).df[:]
    n_genes = len(gene_ontology_term_ids)

    ##
    ## Reduce X
    ##
    big_cube_atts = cube_indexed_dims_no_gene_ontology + cube_non_indexed_dims
    cell_labels, cube_index = make_cube_index(tdb_group, big_cube_atts) # TODO cube_indexed_dims is not an exact match to what was previously passed in, double check that this isnt a breaking change
    n_groups = len(cube_index)

    cube_sum = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_nnz = np.zeros((n_groups, n_genes), dtype=np.uint64)
    cube_min = np.zeros((n_groups, n_genes), dtype=np.float32)
    cube_max = np.zeros((n_groups, n_genes), dtype=np.float32)

    # pass 1 - sum, nnz, min, max
    reduce_X(tdb_group, start_time, cell_labels.cube_idx.values, cube_sum, cube_nnz, cube_min, cube_max)

    return build_in_mem_cube(gene_ontology_term_ids, cube_index, cube_non_indexed_dims, cube_sum, cube_nnz)


def create_empty_cube(uri: str):
    """"
    Create an empty cube with expected schema (dimensions and attributes) at given uri
    """
    if tiledb.array_exists(uri):
        tiledb.remove(uri)

    tiledb.Array.create(
        uri,
        schema,
    )


def create_cube(tdb_group):
    """
    Create queryable cube and write to disk
    """
    uri = f"{tdb_group}/cube"
    start_time = time.time()

    with tiledb.scope_ctx(create_fast_ctx()):
        # Create cube
        create_empty_cube(uri)

        # load data
        dims, vals = load_data_into_cube(tdb_group, uri)

        logger.debug("Saving cube to tiledb")
        with tiledb.open(uri, "w") as cube:
            cube[tuple(dims)] = vals

        logger.debug("Cube created, start consolidation")
        tiledb.consolidate(uri)

        logger.debug("Cube consolidated, start vacuumming")
        tiledb.vacuum(uri)

    logger.debug("Cube creation complete")
    create_cube_sec = time.time() - start_time
    logger.info(f"Big cube: time to create {create_cube_sec}, uri={uri}")
