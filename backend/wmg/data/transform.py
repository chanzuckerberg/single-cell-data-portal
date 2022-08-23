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


def generate_cell_ordering(snapshot_path: str, cell_type_by_tissue: Dict) -> None:
    """
    Use graphviz to map all the cells associated with a tissue to the ontology tree and return their correct order
    """
    # Note: those dependencies are only needed by the WMG pipeline, so we should keep them local
    # so that this file can be imported by tests without breaking.
    from pronto import Ontology
    import pygraphviz as pgv

    onto = Ontology.from_obo_library("cl-basic.obo")

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

        def recurse(node: Set[str], depth=0):

            if node in cells:

                cells.remove(node)
                yield {"id": node, "depth": depth}

                if node != "CL:0000003":
                    depth += 1

            children = [
                (c, positions[c.id]) for c in onto[node].subclasses(with_self=False, distance=1) if c.id in ancestor_ids
            ]
            sorted_children = sorted(children, key=lambda x: x[1][0])
            for child in sorted_children:
                yield from recurse(child[0].id, depth=depth)

        ordered_list = recurse(root)
        return list(ordered_list)

    mapping = {}
    for tissue, cell_df in cell_type_by_tissue.items():
        cells = list(cell_df)
        ordered_cells = compute_ordering(cells, "CL:0000003")
        mapping[tissue] = ordered_cells

    data = []
    for tissue, cells in mapping.items():
        for i, cell in enumerate(cells):
            data.append((tissue, cell["id"], cell["depth"], i))

    df = pd.DataFrame(data, columns=["tissue_ontology_term_id", "cell_type_ontology_term_id", "depth", "order"])
    df.to_json(f"{snapshot_path}/{CELL_TYPE_ORDERINGS_FILENAME}")


def generate_primary_filter_dimensions(snapshot_path: str, corpus_name: str, snapshot_id: int):

    # TODO: remove them from WmgQuery (next 4 following functions)
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

    def order_tissues(ontology_term_ids: Iterable[str]) -> Iterable[str]:
        """
        Order tissues based on appearance in TissueMapper.HIGH_LEVEL_TISSUES
        """
        ontology_term_ids = set(ontology_term_ids)
        ordered_ontology_term_ids = []
        for tissue in TissueMapper.HIGH_LEVEL_TISSUES:
            tissue = TissueMapper.make_id_writable(tissue)
            if tissue in ontology_term_ids:
                ontology_term_ids.remove(tissue)
                ordered_ontology_term_ids.append(tissue)

        if ontology_term_ids:
            ordered_ontology_term_ids = ordered_ontology_term_ids + list(ontology_term_ids)

        return ordered_ontology_term_ids


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
