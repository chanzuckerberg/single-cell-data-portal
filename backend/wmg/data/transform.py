import gc
import logging
import time
from typing import Dict, Union
from anndata._core.views import ArrayView

import numpy
import numpy as np
import tiledb
import pandas as pd
from scipy import sparse
from scipy.sparse import csr_matrix, coo_matrix

from backend.wmg.data.wmg_constants import RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD, INTEGRATED_ARRAY_NAME
from backend.wmg.data.rankit import rankit

from backend.wmg.data.snapshot import CELL_TYPE_ORDERINGS_FILENAME

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def get_cell_types_by_tissue(corpus_group: str) -> Dict:
    """
    Return a list of all associated cell type ontologies for each tissue contained in the
    provided corpus
    """
    with tiledb.open(f"{corpus_group}/obs", "r") as obs:
        tissue_cell_types = (
            obs.query(attrs=[], dims=["tissue_ontology_term_id", "cell_type_ontology_term_id"])
                .df[:]
                .drop_duplicates()
                .sort_values(by="tissue_ontology_term_id")
        )
    unique_tissue_ontology_term_id = tissue_cell_types.tissue_ontology_term_id.unique()
    cell_type_by_tissue = {}
    for x in unique_tissue_ontology_term_id:
        cell_type_by_tissue[x] = tissue_cell_types.loc[
            tissue_cell_types["tissue_ontology_term_id"] == x, "cell_type_ontology_term_id"
        ]

    return cell_type_by_tissue


def generate_cell_ordering(snapshot_path: str, corpus_path: str) -> None:
    """
    Use graphviz to map all the cells associated with a tissue to the ontology tree and return their correct order
    """
    # Note: those dependencies are only needed by the WMG pipeline, so we should keep them local
    # so that this file can be imported by tests without breaking.
    from pronto import Ontology
    import pygraphviz as pgv

    onto = Ontology.from_obo_library("cl-basic.obo")
    cell_type_by_tissue = get_cell_types_by_tissue(corpus_path)
    def compute_ordering(cells, root):
        ancestors = [list(onto[t].superclasses()) for t in cells if t in onto]
        ancestors = [i for s in ancestors for i in s]
        ancestors = set(ancestors)

        G = pgv.AGraph()
        for a in ancestors:
            for s in a.subclasses(with_self=False, distance=1):
                if s in ancestors:
                    G.add_edge(a.id, s.id)

        G.layout(prog="dot")

        positions = {}
        for n in G.iternodes():
            pos = n.attr["pos"].split(",")
            positions[n] = (float(pos[0]), float(pos[1]))

        ancestor_ids = [a.id for a in ancestors]

        def recurse(node):
            if node in cells:
                yield (node)
            children = [
                (c, positions[c.id]) for c in onto[node].subclasses(with_self=False, distance=1) if c.id in ancestor_ids
            ]
            sorted_children = sorted(children, key=lambda x: x[1][0])
            for child in sorted_children:
                yield from recurse(child[0].id)

        ordered_list = list(dict.fromkeys(recurse(root)))
        return ordered_list

    mapping = {}
    for tissue, cell_df in cell_type_by_tissue.items():
        cells = list(cell_df)
        ordered_cells = compute_ordering(cells, "CL:0000003")
        mapping[tissue] = ordered_cells

    data = []
    for tissue, cells in mapping.items():
        for i, cell in enumerate(cells):
            data.append((tissue, cell, i))

    df = pd.DataFrame(data, columns=["tissue_ontology_term_id", "cell_type_ontology_term_id", "order"])
    df.to_json(f"{snapshot_path}/{CELL_TYPE_ORDERINGS_FILENAME}")



# Filter cells
def filter_out_rankits_with_low_expression_counts(
        rankits: csr_matrix, raw_counts_csr: csr_matrix, expect_majority_filtered=True
) -> coo_matrix:
    """
    Keep only rankit values that were computed from expression values above the desired threshold.

    @param expect_majority_filtered: Set this to True if the caller expects the majority of rankit values will be
    filtered, as we can then use an optimal implementation
    """

    rankits_nnz_orig = rankits.nnz
    raw_counts = raw_counts_csr.tocoo(copy=False)
    to_keep_mask = raw_counts.data > RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD

    # Unfortunately, it seems to take longer to decide which strategy to use than it does to just run either one.
    # start = time.time()
    # expect_majority_filtered = to_keep_mask.sum() < rankits_nnz_orig / 2
    # end = time.time()
    # logger.info(f"filter_out_rankits_with_low_expression_counts(): use extract strategy={expect_majority_filtered}; "
    #             f"decision time={end - start}")

    start = time.time()

    if expect_majority_filtered:
        rows_to_keep = raw_counts.row[to_keep_mask]
        cols_to_keep = raw_counts.col[to_keep_mask]
        rankits_to_keep = rankits.data[to_keep_mask]
        rankits_filtered = coo_matrix((rankits_to_keep, (rows_to_keep, cols_to_keep)))
    else:
        to_filter_mask = ~to_keep_mask
        rows_to_filter = raw_counts.row[to_filter_mask]
        cols_to_filter = raw_counts.col[to_filter_mask]
        rankits[rows_to_filter, cols_to_filter] = 0.0
        rankits.eliminate_zeros()
        rankits_filtered = rankits.tocoo()

    end = time.time()

    logger.info(
        f"filter duration={end - start}, "
        f"orig size={rankits_nnz_orig}, "
        f"abs reduction={rankits_nnz_orig - rankits_filtered.nnz}, "
        f"% reduction={(rankits_nnz_orig - rankits_filtered.nnz) / rankits_nnz_orig}"
    )

    return rankits_filtered


# Transform the expression matrix
def transform_expression_raw_counts_to_rankit(
        expression_matrix: Union[np.ndarray, sparse.spmatrix, ArrayView], corpus_path: str, global_var_index: numpy.ndarray, first_obs_idx: int
):
    """
    Apply rankit normalization to raw count expression values and save to the tiledb corpus object
    """
    array_name = f"{corpus_path}/{INTEGRATED_ARRAY_NAME}"
    logger.info(f"saving {array_name}...")
    stride = max(int(np.power(10, np.around(np.log10(1e9 / expression_matrix.shape[1])))), 10_000)
    with tiledb.open(array_name, mode="w") as array:
        for start in range(0, expression_matrix.shape[0], stride):
            end = min(start + stride, expression_matrix.shape[0])
            raw_expression_csr_matrix = sparse.csr_matrix(expression_matrix[start:end, :])

            # Compute RankIt
            rankit_integrated_csr_matrix = rankit(raw_expression_csr_matrix)

            rankit_integrated_coo_matrix = filter_out_rankits_with_low_expression_counts(
                rankit_integrated_csr_matrix, raw_expression_csr_matrix, expect_majority_filtered=True
            )

            global_rows = rankit_integrated_coo_matrix.row + start + first_obs_idx
            global_cols = global_var_index[rankit_integrated_coo_matrix.col]

            rankit_data = rankit_integrated_coo_matrix.data

            assert len(rankit_data) == len(global_rows)
            assert len(rankit_data) == len(global_cols)

            array[global_rows, global_cols] = {"rankit": rankit_data}
            del (
                raw_expression_csr_matrix,
                rankit_integrated_coo_matrix,
                rankit_integrated_csr_matrix,
                global_rows,
                global_cols,
                rankit_data,
            )
            gc.collect()

    logger.debug(f"Saved {array_name}.")


