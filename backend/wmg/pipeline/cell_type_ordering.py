import os
from typing import List, Set

import pandas as pd
import tiledb

from backend.wmg.data.constants import CL_BASIC_OBO_NAME
from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
    CELL_TYPE_ORDERINGS_FILENAME,
)
from backend.wmg.data.utils import get_pinned_ontology_url, to_dict
from backend.wmg.pipeline.constants import (
    CELL_TYPE_ORDERING_CREATED_FLAG,
    EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG,
)
from backend.wmg.pipeline.errors import PipelineStepMissing
from backend.wmg.pipeline.utils import load_pipeline_state, log_func_runtime, write_pipeline_state


@log_func_runtime
def create_cell_type_ordering(corpus_path: str) -> None:
    """
    Writes an ordered "table" of cell types to a json file. The table is written from a pandas data frame with the
    following columns:
        1. tissue_ontology_term_id: str
        2. cell_type_ontology_term_id: str
        3. depth: int  -- 0-based level down in the cell type hierarchy
        4. order: int -- 0-based order per tissue

    :param str snapshot_path: path to write json file
    :param Dict[str, pd.DataFrame] cell_type_by_tissue: all cell types by tissue --
                                                        {"tissue_a": pandas.series with cell types, ...}

    :return None
    """
    pipeline_state = load_pipeline_state(corpus_path)
    if not pipeline_state.get(EXPRESSION_SUMMARY_AND_CELL_COUNTS_CUBE_CREATED_FLAG):
        raise PipelineStepMissing("cell_counts")

    with tiledb.open(os.path.join(corpus_path, CELL_COUNTS_CUBE_NAME)) as cell_counts_cube:
        cell_counts_df = cell_counts_cube.df[:]

    cell_counts_df = (
        cell_counts_df.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"]).first().reset_index()
    )
    cell_type_by_tissue = to_dict(
        cell_counts_df["tissue_ontology_term_id"], cell_counts_df["cell_type_ontology_term_id"].values
    )

    # Generates ordering for ALL cell types
    all_cells = {cell for cell_df in cell_type_by_tissue.values() for cell in cell_df}
    ordered_cells = _cell_type_ordering_compute(all_cells, root="CL:0000003")

    # Create individual data frames per tissue
    ordered_cells_by_tissue = []
    for tissue, target_cells in cell_type_by_tissue.items():
        target_cells = list(target_cells)
        ordered_cells_by_tissue.append(_cell_type_ordering_per_tissue(ordered_cells, target_cells, tissue))

    # Concatenate into single dataframe for writing
    df = pd.concat(ordered_cells_by_tissue).reset_index(drop=True)

    # sort rows for consistency
    df = df.sort_values(by=df.columns.tolist())

    df.to_json(os.path.join(corpus_path, CELL_TYPE_ORDERINGS_FILENAME))

    pipeline_state[CELL_TYPE_ORDERING_CREATED_FLAG] = True
    write_pipeline_state(pipeline_state, corpus_path)


def _cell_type_ordering_compute(cells: Set[str], root: str) -> pd.DataFrame:
    """
    Helper function for `create_cell_type_ordering()`. Orders a set of cell types starting from a root cell type.
    Ordering is according to the CL directed acyclic graph (DAG). Using pygraphviz, a 2-dimensional representation of
    the DAG containing the cell types is built, then an ordered list of cell types is created by traversing the
    DAG using a depth-first approach. Children cell types are represented with a number `depth`, which is 0-based
    and increases as going down through children cell types.

    Any orphan cell types are added at the end of table with depth 0.

    :param Set[str] cells: Set of cell type ontology term ids
    :param str root: Root of the tree, usually CL:0000003

    :return pd.DataFrame: With the following columns:
                             1. cell_type_ontology_term_id: str
                             2. depth: int  -- 0-based level down in the cell type hierarchy.
                                The following nodes are always at depth 0: "CL:0000003", "CL:0000255", "CL:0000548"
    """

    # Note: those dependencies are only needed by the WMG pipeline, so we should keep them local
    # so that this file can be imported by tests without breaking.
    import pygraphviz as pgv
    from pronto import Ontology

    onto = Ontology(get_pinned_ontology_url(CL_BASIC_OBO_NAME))
    ancestors = [list(onto[t].superclasses()) for t in cells if t in onto]
    ancestors = [i for s in ancestors for i in s]
    ancestors = sorted(set(ancestors))

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

    def cell_entity(node, depth):
        return {"id": node, "depth": depth}

    def recurse(node: str, depth=0):
        if node in cells:
            cells.remove(node)
            yield cell_entity(node, depth)

            # Make the root and other high level cell types stay at 0 level. Avoids creating
            # a situation where only the root is at lowest level
            if node not in ["CL:0000003", "CL:0000255", "CL:0000548"]:
                depth += 1

        children = [
            (c, positions[c.id]) for c in onto[node].subclasses(with_self=False, distance=1) if c.id in ancestor_ids
        ]
        sorted_children = sorted(children, key=lambda x: x[1][0])
        for child in sorted_children:
            yield from recurse(child[0].id, depth=depth)

    # Apply recursion to create an ordered list of cells present in set "cells"
    ordered_list = list(recurse(root))

    # If there are any cells left in set "cells", it means that either those cell types
    # don't exist in the ontology or they are above the root ("CL:0000003")
    # Add these "orphan" cells at end of list
    ordered_list.extend([cell_entity(cell, depth=0) for cell in cells])

    # Arranges data into DF
    # Comment these lines out if the desired return object is:
    # [{id:cell_type_1, depth:0}, ... ]
    ordered_df = []
    for cell in ordered_list:
        ordered_df.append([cell["id"], cell["depth"]])
    ordered_df = pd.DataFrame(ordered_df, columns=["cell_type_ontology_term_id", "depth"])

    return ordered_df


def _cell_type_ordering_per_tissue(data_cells: pd.DataFrame, target_cells: List[str], tissue: str) -> pd.DataFrame:
    """
    Helper function for `create_cell_type_ordering()`. Using an ordered cell type table
    obtained by `_cell_type_ordering_compute()`, this function subsets those cell types to only
    contain what's in `target_cells`, and adds a tissue column with values as in `tissue`.

    MOST IMPORTANTLY, this function changes the "depth" column of the ordered cell type table based
    on the `target_cells`

    :param pd.DataFrame data_cells: table of ordered cell types from `_cell_type_ordering_compute()`
    :param target_cells List[str]: cell types of interest.
    :param str tissue: Tissue label, it will be added as column to final table

    :return pd.Dataframe: with the following columns:
                             1. tissue_ontology_term_id: str -- as set in `tissue` parameter
                             2. cell_type_ontology_term_id: str -- only with cells from `target_cells` parameter
                             3. depth: int -- fixed depths
                             4. order: int -- 0:row_numbers
    """

    tissue_cells = data_cells.copy()
    tissue_cells["is_in_tissue"] = [cell in target_cells for cell in tissue_cells["cell_type_ontology_term_id"]]

    depth_col = tissue_cells.columns.get_loc("depth")

    # Place holder column for row that belong to tissue
    is_in_tissue_col = tissue_cells.columns.get_loc("is_in_tissue")

    # Re arranged depths based on the rows that belong to tissue
    for i in range(len(tissue_cells)):
        if not tissue_cells.iloc[i, is_in_tissue_col]:
            original_depth = tissue_cells.iloc[i, depth_col]
            for j in range(i + 1, len(tissue_cells)):
                if original_depth < tissue_cells.iloc[j, depth_col]:
                    tissue_cells.iloc[j, depth_col] -= 1
                else:
                    break

    # Remove placeholder column
    tissue_cells = tissue_cells[tissue_cells["is_in_tissue"]]
    del tissue_cells["is_in_tissue"]

    # Insert first column as tissue
    tissue_cells.insert(0, "tissue_ontology_term_id", tissue)
    # Insert last column as order
    tissue_cells["order"] = range(len(tissue_cells))

    return tissue_cells
