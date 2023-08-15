import json
from typing import Dict, Iterable, List, Set

import pandas as pd
import tiledb

from backend.wmg.data.ontology_labels import gene_term_label, ontology_term_label
from backend.wmg.data.snapshot import (
    CELL_TYPE_ORDERINGS_FILENAME,
    EXPRESSION_SUMMARY_CUBE_NAME,
    PRIMARY_FILTER_DIMENSIONS_FILENAME,
)
from backend.wmg.data.tissue_mapper import TissueMapper
from backend.wmg.data.utils import get_pinned_ontology_url, log_func_runtime, to_dict


@log_func_runtime
def cell_type_ordering_create_file(snapshot_path: str) -> None:
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

    cell_counts_cube = tiledb.open(f"{snapshot_path}/cell_counts")
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

    df.to_json(f"{snapshot_path}/{CELL_TYPE_ORDERINGS_FILENAME}")


def _cell_type_ordering_compute(cells: Set[str], root: str) -> pd.DataFrame:
    """
    Helper function for `cell_type_ordering_create_file()`. Orders a set of cell types starting from a root cell type.
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

    onto = Ontology(get_pinned_ontology_url("cl-basic.obo"))
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
    Helper function for `cell_type_ordering_create_file()`. Using an ordered cell type table
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


@log_func_runtime
def generate_primary_filter_dimensions(snapshot_path: str, snapshot_id: int):
    def list_primary_filter_dimension_term_ids(cube, primary_dim_name: str):
        return cube.query(attrs=[], dims=[primary_dim_name]).df[:].groupby([primary_dim_name]).first().index.tolist()

    def list_grouped_primary_filter_dimensions_term_ids(
        cube, primary_dim_name: str, group_by_dim: str
    ) -> Dict[str, List[str]]:
        return (
            cube.query(attrs=[], dims=[primary_dim_name, group_by_dim])
            .df[:]
            .drop_duplicates()
            .groupby(group_by_dim)
            .agg(list)
            .to_dict()[primary_dim_name]
        )

    def build_gene_id_label_mapping(gene_ontology_term_ids: List[str]) -> List[dict]:
        return [
            {gene_ontology_term_id: gene_term_label(gene_ontology_term_id)}
            for gene_ontology_term_id in gene_ontology_term_ids
        ]

    def build_ontology_term_id_label_mapping(ontology_term_ids: Iterable[str]) -> List[dict]:
        return [{ontology_term_id: ontology_term_label(ontology_term_id)} for ontology_term_id in ontology_term_ids]

    with tiledb.open(f"{snapshot_path}/{EXPRESSION_SUMMARY_CUBE_NAME}") as cube:

        # gene terms are grouped by organism, and represented as a nested lists in dict, keyed by organism
        organism_gene_ids: dict[str, List[str]] = list_grouped_primary_filter_dimensions_term_ids(
            cube, "gene_ontology_term_id", group_by_dim="organism_ontology_term_id"
        )
        organism_gene_terms = {
            organism_term_id: build_gene_id_label_mapping(gene_term_ids)
            for organism_term_id, gene_term_ids in organism_gene_ids.items()
        }

        # tissue terms are grouped by organism, and represented as a nested lists in dict, keyed by organism
        organism_tissue_ids: dict[str, List[str]] = list_grouped_primary_filter_dimensions_term_ids(
            cube, "tissue_ontology_term_id", group_by_dim="organism_ontology_term_id"
        )
        organism_tissue_terms = {
            organism_term_id: build_ontology_term_id_label_mapping(order_tissues(tissue_term_ids))
            for organism_term_id, tissue_term_ids in organism_tissue_ids.items()
        }

        result = dict(
            snapshot_id=str(snapshot_id),
            organism_terms=build_ontology_term_id_label_mapping(
                list_primary_filter_dimension_term_ids(cube, "organism_ontology_term_id")
            ),
            tissue_terms=organism_tissue_terms,
            gene_terms=organism_gene_terms,
        )

        with open(f"{snapshot_path}/{PRIMARY_FILTER_DIMENSIONS_FILENAME}", "w") as f:
            json.dump(result, f)


def order_tissues(ontology_term_ids: Iterable[str]) -> Iterable[str]:
    """
    Order tissues based on appearance in TissueMapper.HIGH_LEVEL_TISSUES. This will maintain the priority set in
    that class which is intended to keep most relevant tissues on top and tissues that are related to be placed
    sequentially
    """
    ontology_term_ids = set(ontology_term_ids)
    ordered_ontology_term_ids = []
    for tissue in TissueMapper.HIGH_LEVEL_TISSUES:
        tissue = TissueMapper.reformat_ontology_term_id(tissue, to_writable=True)
        if tissue in ontology_term_ids:
            ontology_term_ids.remove(tissue)
            ordered_ontology_term_ids.append(tissue)

    if ontology_term_ids:
        ordered_ontology_term_ids = ordered_ontology_term_ids + list(ontology_term_ids)

    return ordered_ontology_term_ids
