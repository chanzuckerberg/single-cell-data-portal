import logging

import pandas as pd
from pronto import Ontology

from backend.common.utils.rollup import rollup_across_cell_type_descendants
from backend.wmg.data.constants import CL_BASIC_PERMANENT_URL_PRONTO

logger = logging.getLogger(__name__)


"""
This module contains the OntologyTreeBuilder class which is used to build a nested dictionary representation of the cell ontology tree.
The tree is populated with cell type counts from an input dataframe.
"""


class OntologyTreeBuilder:
    def __init__(self, cell_counts_df, root_node="CL:0000548"):
        """
        OntologyTreeBuilder is a class that builds a nested dictionary representation of the cell ontology tree
        and populates the nodes with cell type counts from the input cell counts dataframe.

        Arguments
        ---------
        cell_counts_df: pandas.DataFrame
            A dataframe with required columns "tissue_ontology_term_id", "cell_type_ontology_term_id", and "n_cells".
            This dataframe is typically the cell counts cube from the WMG snapshot.

        root_node: str, optional, default="CL:0000548" (animal cell)
            The root node of the ontology tree. This is the node from which the ontology tree is traversed.
        """

        logger.info("Loading CL ontology...")
        self.ontology = Ontology(CL_BASIC_PERMANENT_URL_PRONTO)

        logger.info("Initializing cell type data structures from the input cell counts dataframe...")
        self.id_to_name, self.all_cell_type_ids = self._initialize_id_to_name()
        self.cell_counts_df, self.cell_counts_df_rollup = self._initialize_cell_counts_df_rollup(cell_counts_df)
        self.all_cell_type_ids_in_corpus = self.cell_counts_df_rollup.index.values[
            self.cell_counts_df_rollup.values > 0
        ]
        logger.info("Initializing ontology tree data structures by traversing CL ontology...")
        self.ontology_graph, self.traverse_node_counter, self.all_unique_nodes = self._traverse_ontology_with_counting(
            self.ontology[root_node]
        )
        self._delete_unknown_terms_from_ontology_graph(self.ontology_graph)
        self.start_node = f"{root_node}__0"

    def _initialize_cell_counts_df_rollup(self, cell_counts_df):
        """
        This function prepares the cell counts dataframe and its rollup for the ontology tree traversal.

        Arguments
        ---------
        cell_counts_df: pandas.DataFrame
            This dataframe is typically the cell counts cube from the WMG snapshot.

        Returns
        -------
        cell_counts_df: pandas.Series
            A series with index "cell_type_ontology_term_id" and values "n_cells".
            This series is the cell counts dataframe with the following modification:
            1. The cell type ontology term ids that are not in the input cell counts dataframe are added with 0 counts.

        cell_counts_df_rollup: pandas.Series
            A series with index "cell_type_ontology_term_id" and values "n_cells".
            "n_cells" contains the sum of "n_cells" for the cell type and its descendants.
        """

        cell_counts_df = (
            cell_counts_df.groupby("cell_type_ontology_term_id").sum(numeric_only=True)[["n_cells"]].reset_index()
        )

        to_attach = pd.DataFrame()
        to_attach["cell_type_ontology_term_id"] = [
            i for i in self.all_cell_type_ids if i not in cell_counts_df["cell_type_ontology_term_id"].values
        ]
        to_attach["n_cells"] = 0
        cell_counts_df = pd.concat([cell_counts_df, to_attach], axis=0)
        cell_counts_df_rollup = rollup_across_cell_type_descendants(cell_counts_df)
        cell_counts_df = cell_counts_df.set_index("cell_type_ontology_term_id")["n_cells"]
        cell_counts_df_rollup = cell_counts_df_rollup.set_index("cell_type_ontology_term_id")["n_cells"]
        return cell_counts_df, cell_counts_df_rollup

    def _delete_unknown_terms_from_ontology_graph(self, node):
        """
        Deletes nodes that have no name in the CL ontology from the ontology graph.

        Arguments
        ---------
        node: dict
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        """
        new_children = []
        for child in node.get("children", []):
            unknown = child["name"].startswith("CL:")
            if not unknown:
                new_children.append(child)
        if len(new_children) > 0:
            node["children"] = new_children
        elif "children" in node:
            del node["children"]

        for child in node.get("children", []):
            self._delete_unknown_terms_from_ontology_graph(child)

    def _initialize_id_to_name(self):
        """
        Initializes a dictionary that maps cell type ontology term ids to their names.

        Returns
        -------
        id_to_name: dict
            A dictionary that maps cell type ontology term ids to their names.
        all_cell_type_ids: list
            A list of all cell type ontology term ids.
        """
        id_to_name = {}
        for term in self.ontology.terms():
            if term.id.startswith("CL:"):
                id_to_name[term.id] = term.name
        all_cell_type_ids = list(id_to_name.keys())
        return id_to_name, all_cell_type_ids

    def _traverse_ontology_with_counting(self, node, traverse_node_counter=None, all_unique_nodes=None):
        """
        This function traverses the cell ontology and builds a nested dictionary representation of the tree.
        It also counts the number of times a node has been visited and adds a suffix to the node id to make it unique.
        The suffix is the number of times the node had been visited prior.

        Arguments
        ---------
        node: pronto.Term
            A node of the cell ontology. This should be the root node of the ontology at the top-level.
        traverse_node_counter: dict, optional, default=None
            A dictionary that maps cell type ontology term ids to the number of times they have been visited.
        all_unique_nodes: set, optional, default=None
            A set of all unique cell type ontology term ids that have been visited.

        Returns
        -------
        ontology_graph: dict
            A nested dictionary representation of the cell ontology tree.
        """
        if traverse_node_counter is None:
            traverse_node_counter = {}
        if all_unique_nodes is None:
            all_unique_nodes = set()

        node_count = traverse_node_counter.get(node.id, 0)
        traverse_node_counter[node.id] = node_count + 1
        all_unique_nodes.add(node.id + "__" + str(node_count))

        subclasses = list(node.subclasses(with_self=False, distance=1))

        if len(subclasses) == 0:
            return (
                {
                    "id": node.id + "__" + str(node_count),
                    "name": self.id_to_name[node.id] if node.id in self.id_to_name else node.id,
                    "n_cells_rollup": int(
                        self.cell_counts_df_rollup[node.id] if node.id in self.cell_counts_df_rollup else 0
                    ),
                    "n_cells": int(self.cell_counts_df[node.id] if node.id in self.cell_counts_df else 0),
                },
                None,
                None,
            )

        children = []
        for child in subclasses:
            subtree, _, _ = self._traverse_ontology_with_counting(child, traverse_node_counter, all_unique_nodes)
            children.append(subtree)

        return (
            {
                "id": node.id + "__" + str(node_count),
                "name": self.id_to_name[node.id] if node.id in self.id_to_name else node.id,
                "n_cells_rollup": int(
                    self.cell_counts_df_rollup[node.id] if node.id in self.cell_counts_df_rollup else 0
                ),
                "n_cells": int(self.cell_counts_df[node.id] if node.id in self.cell_counts_df else 0),
                "children": children,
            },
            traverse_node_counter,
            all_unique_nodes,
        )
