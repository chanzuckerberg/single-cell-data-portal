import json
import logging
from dataclasses import dataclass
from typing import Any, Dict, Optional

import pandas as pd
from pronto import Ontology, Term

from backend.common.utils.rollup import rollup_across_cell_type_descendants
from backend.wmg.data.constants import CL_BASIC_PERMANENT_URL_PRONTO

logger = logging.getLogger(__name__)


"""
This module contains the OntologyTreeBuilder class which is used to build a nested dictionary representation of the cell ontology tree.
The tree is populated with cell type counts from an input dataframe.
"""


@dataclass
class TraverseOntologyResult:
    subtree: Dict[str, Any]
    traverse_node_counter: Dict[str, int]
    all_unique_nodes: set[str]


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

        logger.info(f"Loading CL ontology from root node {root_node}...")
        self.ontology = Ontology(CL_BASIC_PERMANENT_URL_PRONTO)

        logger.info("Initializing cell type data structures from the input cell counts dataframe...")
        self.id_to_name, self.all_cell_type_ids = self._initialize_id_to_name()
        self.cell_counts_df, self.cell_counts_df_rollup = self._initialize_cell_counts_df_rollup(cell_counts_df)
        self.all_cell_type_ids_in_corpus = self.cell_counts_df_rollup.index.values[
            self.cell_counts_df_rollup.values > 0
        ]
        logger.info("Initializing ontology tree data structures by traversing CL ontology...")
        traverse_ontology_result = self._traverse_ontology_with_counting(self.ontology[root_node])
        self.ontology_graph = traverse_ontology_result.subtree
        self.traverse_node_counter = traverse_ontology_result.traverse_node_counter
        self.all_unique_nodes = traverse_ontology_result.all_unique_nodes

        self._delete_unknown_terms_from_ontology_graph(self.ontology_graph)

        self.all_children, self.all_parents = self._initialize_children_and_parents_per_node()
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

        # to_attach is a DataFrame that will contain cell type ontology term ids that are not present in the input cell counts dataframe. These will be added with 0 counts.
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

    def _traverse_ontology_with_counting(
        self,
        node: Term,
        traverse_node_counter: Optional[Dict[str, int]] = None,
        all_unique_nodes: Optional[set[str]] = None,
    ) -> TraverseOntologyResult:
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
            return TraverseOntologyResult(
                subtree={
                    "id": node.id + "__" + str(node_count),
                    "name": self.id_to_name[node.id] if node.id in self.id_to_name else node.id,
                    "n_cells_rollup": int(
                        self.cell_counts_df_rollup[node.id] if node.id in self.cell_counts_df_rollup else 0
                    ),
                    "n_cells": int(self.cell_counts_df[node.id] if node.id in self.cell_counts_df else 0),
                },
                traverse_node_counter=None,
                all_unique_nodes=None,
            )

        children = []
        for child in subclasses:
            traverse_ontology_result = self._traverse_ontology_with_counting(
                child, traverse_node_counter, all_unique_nodes
            )
            children.append(traverse_ontology_result.subtree)

        return TraverseOntologyResult(
            subtree={
                "id": node.id + "__" + str(node_count),
                "name": self.id_to_name[node.id] if node.id in self.id_to_name else node.id,
                "n_cells_rollup": int(
                    self.cell_counts_df_rollup[node.id] if node.id in self.cell_counts_df_rollup else 0
                ),
                "n_cells": int(self.cell_counts_df[node.id] if node.id in self.cell_counts_df else 0),
                "children": children,
            },
            traverse_node_counter=traverse_node_counter,
            all_unique_nodes=all_unique_nodes,
        )

    ### Get the ontology tree state
    def get_ontology_tree_state_per_celltype(self):
        """
        This function gets the ontology tree state per cell type. The ontology tree state is a mask that determines
        which nodes are expanded by default and which nodes are not shown when expanded in the ontology.

        The following rules determine which cell types are expanded by default and not shown when expanded:
         - All instances of the target cell type are shown.
         - All children of the first instance of the target cell type are shown.
         - All grandchildren of the first instance of the target cell type are shown.
         - All siblings of the first instance of the target cell type are shown.
         - All nodes in the paths to any instance of the target cell type are shown.
        Nodes that are not shown when expanded will be collapsed into a dummy node that contains the ids of the
        children that are not shown. Nodes will be collapsed into dummy nodes if they are siblings of nodes that are
        in a path to any instances of the target cell type.

        Returns
        -------
        all_states_per_cell_type: dict
            A dictionary that maps cell type ontology term ids to their ontology tree states.
            The ontology tree state is a dictionary with keys
             - "isExpandedNodes": A list of cell type ontology term ids that are expanded by default.
             - "notShownWhenExpandedNodes": A dictionary that maps cell type ontology term ids to their hidden children
        """
        all_states_per_cell_type = {}
        for i, end_node in enumerate(self.all_cell_type_ids_in_corpus):
            if end_node in self.traverse_node_counter:
                logger.info(
                    "Getting ontology tree state for cell type %s (%s/%s)",
                    end_node,
                    i + 1,
                    len(self.all_cell_type_ids_in_corpus),
                )

                all_paths = []
                siblings = []
                for i in range(self.traverse_node_counter[end_node]):
                    end_node_i = end_node + "__" + str(i)
                    path = self._depth_first_search_pathfinder(end_node_i)
                    path = path if path else [end_node_i]
                    all_paths.append(path)

                    siblings.append(
                        sum([self.all_children.get(parent, []) for parent in self.all_parents.get(end_node_i, [])], [])
                    )

                visited_nodes_in_paths = list(set(sum(all_paths, [])))  # in a path to target

                children = self.all_children.get(end_node + "__0", [])  # children
                grandchildren = sum([self.all_children.get(child, []) for child in children], [])  # grandchildren
                siblings = list(set(sum(siblings, [])))  # siblings
                valid_nodes = list(set(visited_nodes_in_paths + children + grandchildren + siblings))

                ontology_graph_copy = self._get_deepcopy_of_ontology_graph()
                self._truncate_graph_first_pass(ontology_graph_copy, valid_nodes)
                self._truncate_graph_second_pass(ontology_graph_copy, visited_nodes_in_paths)

                isExpandedNodes = list(set(_getExpandedData(ontology_graph_copy)))
                notShownWhenExpandedNodes = _getShownData(ontology_graph_copy)

                notShownWhenExpanded = {}
                for i in notShownWhenExpandedNodes:
                    notShownWhenExpanded.update(i)
                all_states_per_cell_type[end_node] = {
                    "isExpandedNodes": isExpandedNodes,
                    "notShownWhenExpandedNodes": notShownWhenExpanded,
                }

        return all_states_per_cell_type

    def _depth_first_search_pathfinder(self, path_end_node, node=None, path=None):
        """
        This function finds a path backwards from the end node to the start node using depth-first search.

        Arguments
        ---------
        path_end_node: str
            The end node of the path.
        node: str, optional, default=None
            The current node in the depth-first search.
        path: list, optional, default=None
            The current path in the depth-first search.

        Returns
        -------
        path: list
            The path from the end node to the start node.
        """
        if path is None:
            path = [path_end_node]
            node = path_end_node

        if node == self.start_node:
            return path

        for parent in self.all_parents.get(node, []):
            full_path = self._depth_first_search_pathfinder(path_end_node, node=parent, path=path + [parent])
            if full_path:
                return full_path

    def _truncate_graph_first_pass(self, graph, valid_nodes):
        """
        This function truncates the ontology graph in-place by removing nodes that are not in the valid nodes list.
        It also adds a dummy node to the graph if the graph has children that are not in the valid nodes list.
        The dummy node contains the ids of the invalid children.

        Arguments
        ---------
        graph: dict
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        valid_nodes: list
            A list of valid cell type ontology term ids.
        """
        if graph["id"] not in valid_nodes:
            return False

        children = graph.get("children", [])
        valid_children = []
        append_dummy = False

        invalid_children_ids = []
        for child in children:
            is_valid = self._truncate_graph_first_pass(child, valid_nodes)
            if is_valid:
                valid_children.append(child)
            elif child["id"] != "":
                invalid_children_ids.append(child["id"])
                append_dummy = True

        if append_dummy and len(valid_children) > 0:
            valid_children.append(
                {
                    "id": "",
                    "name": "",
                    "n_cells_rollup": 0,
                    "n_cells": 0,
                    "invalid_children_ids": invalid_children_ids,
                    "parent": graph["id"],
                }
            )
        if len(valid_children) > 0:
            graph["children"] = valid_children
        else:
            if "children" in graph:
                del graph["children"]

        return True

    def _truncate_graph_second_pass(self, graph, visited_nodes_in_paths, nodesWithChildrenFound=None):
        """
        We do not want to show cell types multiple times in the ontology tree unless they are nodes that are in one of the
        paths to the target cell type. In this case, show the node in the path and collapse all its siblings into the
        dummy node.

        Arguments
        ---------
        graph: dict
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        visited_nodes_in_paths: list
            A list of cell type ontology term ids that are in one of the paths to the target cell type.
            The target cell type may have multiple instances in the ontology graph (each instance with a unique suffix).
        nodesWithChildrenFound: set, optional, default=None
            A set of cell type ontology term ids that have children. This is used to remove children of nodes that
            have already been encountered during the traversal.
        """
        if nodesWithChildrenFound is None:
            nodesWithChildrenFound = set()
        if graph["id"].split("__")[0] in nodesWithChildrenFound:
            if "children" in graph:
                children = graph["children"]
                new_children = []
                invalid_children_ids = []
                for child in children:
                    if child["id"] in visited_nodes_in_paths:
                        new_children.append(child)
                    elif child["id"] != "":
                        invalid_children_ids.append(child["id"])

                if len(children) > len(new_children) and len(new_children) > 0:
                    # append dummy
                    new_children.append(
                        {
                            "id": "",
                            "name": "",
                            "n_cells_rollup": 0,
                            "n_cells": 0,
                            "invalid_children_ids": invalid_children_ids,
                            "parent": graph["id"],
                        }
                    )
                if len(new_children) > 0:
                    graph["children"] = new_children
                else:
                    del graph["children"]
        elif "children" in graph:
            nodesWithChildrenFound.add(graph["id"].split("__")[0])

        children = graph.get("children", [])
        for child in children:
            if child["id"] != "":
                self._truncate_graph_second_pass(child, visited_nodes_in_paths, nodesWithChildrenFound)

    def _initialize_children_and_parents_per_node(self):
        """
        Initializes dictionaries that map cell type ontology term ids to their children and parents.
        While the cell ontology contains multiple instances of the same node, the ontology graph
        contains only one instance of each node as each cell type has been suffixed with a unique integer.
        Therefore, each node only has one parent and a unique set of children.

        Returns
        -------
        all_children: dict
            A dictionary that maps cell type ontology term ids to their children.
        all_parents: dict
            A dictionary that maps cell type ontology term ids to their parents.
        """
        all_children = self._build_children_per_node(self.ontology_graph)
        all_parents = self._build_parents_per_node(self.ontology_graph)
        return all_children, all_parents

    def _build_children_per_node(self, node, all_children=None):
        """
        Builds a dictionary that maps cell type ontology term ids to their children.

        Arguments
        ---------
        node: dict
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        all_children: dict, optional, default=None
            A dictionary that maps cell type ontology term ids to their children.

        Returns
        -------
        all_children: dict
            A dictionary that maps cell type ontology term ids to their children.
        """
        if all_children is None:
            all_children = {}
        children = node.get("children", [])
        ids = [] if len(children) == 0 else [child["id"] for child in children]

        all_children[node["id"]] = ids
        for child in children:
            self._build_children_per_node(child, all_children)
        return all_children

    def _build_parents_per_node(self, node, all_parents=None):
        """
        Builds a dictionary that maps cell type ontology term ids to their parents.

        Arguments
        ---------
        node: dict
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        all_parents: dict, optional, default=None
            A dictionary that maps cell type ontology term ids to their parents.

        Returns
        -------
        all_parents: dict
            A dictionary that maps cell type ontology term ids to their parents.
        """
        if all_parents is None:
            all_parents = {}
        children = node.get("children", [])

        for child in children:
            all_parents[child["id"]] = [node["id"]]
            self._build_parents_per_node(child, all_parents)
        return all_parents

    def _get_deepcopy_of_ontology_graph(self):
        """
        Returns a deepcopy of the ontology graph.
        """
        return json.loads(json.dumps(self.ontology_graph))


def _getExpandedData(ontology_graph, isExpandedNodes=None):
    """
    This function gets the cell type ontology term ids of the nodes that are expanded by default.

    Arguments
    ---------
    ontology_graph: dict
        A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.

    Returns
    -------
    isExpandedNodes: list
        A list of cell type ontology term ids that are expanded by default.
    """
    if isExpandedNodes is None:
        isExpandedNodes = []

    if "children" in ontology_graph:
        isExpandedNodes.append(ontology_graph["id"])
        for child in ontology_graph["children"]:
            _getExpandedData(child, isExpandedNodes)

    return isExpandedNodes


def _getShownData(graph, notShownWhenExpandedNodes=None):
    """
    This function gets the children nodes that are not shown when expanded for each node in the ontology graph.

    Arguments
    ---------
    graph: dict
        A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.

    Returns
    -------
    notShownWhenExpandedNodes: list
        A list of dictionaries that map cell type ontology term ids to their hidden children.
    """
    if notShownWhenExpandedNodes is None:
        notShownWhenExpandedNodes = []

    if "children" in graph:
        for child in graph["children"]:
            if child["id"] == "":
                if len(child["invalid_children_ids"]) > 0:
                    notShownWhenExpandedNodes.append({child["parent"]: list(set(child["invalid_children_ids"]))})
            else:
                _getShownData(child, notShownWhenExpandedNodes)
    return notShownWhenExpandedNodes
