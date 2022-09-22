import tiledb
import pandas as pd

from backend.wmg.data.ontology_labels import ontology_term_label, gene_term_label
from backend.wmg.data.tissue_mapper import TissueMapper
from typing import Dict, List, Iterable, Set
import json

from backend.wmg.data.snapshot import (
    CELL_TYPE_ORDERINGS_FILENAME,
    EXPRESSION_SUMMARY_CUBE_NAME,
    PRIMARY_FILTER_DIMENSIONS_FILENAME,
)
from backend.wmg.data.utils import log_func_runtime


@log_func_runtime
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


@log_func_runtime
def generate_cell_ordering(snapshot_path: str, cell_type_by_tissue: Dict) -> None:
    """
    Use graphviz to map all the cells associated with a tissue to the ontology tree and return their correct order
    """
    # Note: those dependencies are only needed by the WMG pipeline, so we should keep them local
    # so that this file can be imported by tests without breaking.
    from pronto import Ontology
    import pygraphviz as pgv

    onto = Ontology.from_obo_library("cl-basic.obo")

    def compute_ordering(cells: Set[str], root: str):
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

                # Make the root stay in the same level. Avoids creating
                # a situation where only the root is at lowest level
                if node != root:
                    depth += 1

            children = [
                (c, positions[c.id]) for c in onto[node].subclasses(with_self=False, distance=1) if c.id in ancestor_ids
            ]
            sorted_children = sorted(children, key=lambda x: x[1][0])
            for child in sorted_children:
                yield from recurse(child[0].id, depth=depth)

        # Apply recursion to create ordered list of cells present in set "cells"
        ordered_list = list(recurse(root))

        # If there are any cells left in set "cells", it means that either those cell types
        # don't exist in the ontology or they are above the root ("CL:0000003")
        # Add these "orphan" cells at end of list
        while cells:
            ordered_list.append(cell_entity(cells.pop(), depth=0))

        # Arranges data into DF
        # Comment these lines out if the desired return object is:
        # [{id:cell_type_1, depth:0}, ... ]
        ordered_df = []
        for cell in ordered_list:
            ordered_df.append([cell["id"], cell["depth"]])
        ordered_df = pd.DataFrame(ordered_df, columns=["cell_type_ontology_term_id", "depth"])
        return ordered_df

    def create_tissue_cell_type_df(data_cells: pd.DataFrame, target_cells: list, tissue):

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

    # Generates ordering for ALL cell types
    all_cells = {cell for cell_df in cell_type_by_tissue.values() for cell in cell_df}
    ordered_cells = compute_ordering(all_cells, "CL:0000003")

    # Create individual data frames per tissue
    ordered_cells_by_tissue = []
    for tissue, target_cells in cell_type_by_tissue.items():
        target_cells = list(target_cells)
        ordered_cells_by_tissue.append(create_tissue_cell_type_df(ordered_cells, target_cells, tissue))

    df = pd.concat(ordered_cells_by_tissue).reset_index(drop=True)

    df.to_json(f"{snapshot_path}/{CELL_TYPE_ORDERINGS_FILENAME}")




@log_func_runtime
def generate_primary_filter_dimensions(snapshot_path: str, corpus_name: str, snapshot_id: int):
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

    with tiledb.open(f"{snapshot_path}/{corpus_name}/{EXPRESSION_SUMMARY_CUBE_NAME}") as cube:

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
