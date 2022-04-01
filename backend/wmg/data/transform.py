from typing import Dict
import tiledb
import pandas as pd

from backend.wmg.data.snapshot import CELL_TYPE_ORDERINGS_FILENAME
from pronto import Ontology
import pygraphviz as pgv


def get_cells_by_tissue_type(corpus_group: str) -> Dict:
    """
    Return a list of all associated cell type ontologies for each tissue contained in the
    provided corpus
    """
    with tiledb.open(f"{corpus_group}/obs", "r") as obs:
        cell_tissue_types = (
            obs.query(attrs=[], dims=["tissue_ontology_term_id", "cell_type_ontology_term_id"])
            .df[:]
            .drop_duplicates()
            .sort_values(by="tissue_ontology_term_id")
        )
    unique_tissue_ontology_term_id = cell_tissue_types.tissue_ontology_term_id.unique()
    cell_type_by_tissue = {}
    for x in unique_tissue_ontology_term_id:
        cell_type_by_tissue[x] = cell_tissue_types.loc[
            cell_tissue_types["tissue_ontology_term_id"] == x, "cell_type_ontology_term_id"
        ]

    return cell_type_by_tissue


def generate_cell_ordering(cell_type_by_tissue: Dict) -> None:
    """
    Use graphviz to map all the cells assoicated with a tissue to the ontology tree and return their correct order
    """
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
    df.to_json(CELL_TYPE_ORDERINGS_FILENAME)
