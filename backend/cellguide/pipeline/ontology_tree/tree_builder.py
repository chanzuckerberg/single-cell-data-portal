import copy
import logging
import warnings
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd
from dask import compute, delayed
from dask.diagnostics import ProgressBar
from pronto import Ontology, Term

from backend.cellguide.pipeline.constants import CELLGUIDE_PIPELINE_NUM_CPUS, UBERON_BASIC_PERMANENT_URL_PRONTO
from backend.cellguide.pipeline.ontology_tree.types import OntologyTree, OntologyTreeState
from backend.common.utils.rollup import rollup_across_cell_type_descendants
from backend.wmg.data.constants import CL_BASIC_OBO_NAME
from backend.wmg.data.utils import get_pinned_ontology_url

logger = logging.getLogger(__name__)


"""
This module contains the OntologyTreeBuilder class which is used to build a nested dictionary representation of the cell ontology tree.
The tree is populated with cell type counts from an input dataframe. The class also provides functions to get the ontology tree state
per cell type and per tissue. The ontology tree state is a mask that determines which nodes are expanded by default and which nodes 
are shown when expanded in the ontology.

The module also defines a list of tissues (TISSUES_IMMUNE_CELL_WHITELIST) in which we want to show immune cell types and a constant 
HEMATOPOIETIC_CELL_TYPE_ID representing the immune cell type ancestor to be filtered out from non-whitelisted tissues.
"""

# these nodes will be traversed first compared to their siblings
# this is to ensure that animal cell descendants are encountered first
# such that {cell_type}__0 is always a descendant of animal cell if
# {cell_type} is a descendant of animal cell.
TRAVERSAL_PRIORITY_NODES = [
    "CL:0000000",  # cell
    "CL:0000003",  # native cell,
    "CL:0000255",  # eukaryotic cell
    "CL:0000548",  # animal cell
]

# tissues in which we want to show immune cell types
TISSUES_IMMUNE_CELL_WHITELIST = [
    "UBERON:0002193",  # hemolymphoid system
    "UBERON:0002390",  # hematopoietic system
    "UBERON:0000178",  # blood
    "UBERON:0005057",  # immune organ
]
# immune cell type ancestor to be filtered out from non-whitelisted tissues
HEMATOPOIETIC_CELL_TYPE_ID = "CL:0000988"


@dataclass
class TraverseOntologyResult:
    subtree: OntologyTree
    traverse_node_counter: Dict[str, int]
    all_unique_nodes: set[str]


class OntologyTreeBuilder:
    def __init__(self, cell_counts_df, root_node="CL:0000000"):
        """
        OntologyTreeBuilder is a class that builds a nested dictionary representation of the cell ontology tree
        and populates the nodes with cell type counts from the input cell counts dataframe. It also provides
        functions to get the ontology tree state per cell type and per tissue. The ontology tree state is a
        mask that determines which nodes are expanded by default and which nodes are shown when expanded in the
        ontology.

        Arguments
        ---------
        cell_counts_df: pandas.DataFrame
            A dataframe with required columns "tissue_ontology_term_id", "cell_type_ontology_term_id", and "n_cells".
            This dataframe is typically the cell counts cube from the WMG snapshot.

        root_node: str, optional, default="CL:0000548" (animal cell)
            The root node of the ontology tree. This is the node from which the ontology tree is traversed.
        """

        logger.info(f"Loading CL ontology from root node {root_node}...")
        self.ontology = Ontology(get_pinned_ontology_url(CL_BASIC_OBO_NAME))

        logger.info("Loading UBERON ontology...")
        with warnings.catch_warnings():
            # loading uberon ontology has some warnings that we don't care about
            warnings.simplefilter("ignore")
            self.uberon_ontology = Ontology(UBERON_BASIC_PERMANENT_URL_PRONTO)

        logger.info("Initializing tissue data structures from the input cell counts dataframe...")
        self.tissue_counts_df = cell_counts_df.groupby("tissue_ontology_term_id").sum(numeric_only=True)["n_cells"]
        self.tissue_celltypes_df = (
            cell_counts_df.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"])
            .sum(numeric_only=True)
            .reset_index()
        )
        self.uberon_by_celltype = _to_dict(
            cell_counts_df["tissue_ontology_term_id"], cell_counts_df["cell_type_ontology_term_id"]
        )

        logger.info("Initializing cell type data structures from the input cell counts dataframe...")
        self.id_to_name, self.all_cell_type_ids = self._initialize_id_to_name()
        self.cell_counts_df, self.cell_counts_df_rollup = self._initialize_cell_counts_df_rollup(cell_counts_df)
        self.all_cell_type_ids_in_corpus = self.cell_counts_df_rollup.index.values[
            self.cell_counts_df_rollup.values > 0
        ]

        # used for gpt pipeline
        # TODO: refactor this module to use the below mapping instead of self.id_to_name
        self.all_cell_type_ids_to_labels_in_corpus = dict(
            zip(self.all_cell_type_ids_in_corpus, [self.id_to_name.get(c, c) for c in self.all_cell_type_ids_in_corpus])
        )

        self.all_tissue_ids_in_corpus = list(self.uberon_by_celltype.keys())
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

    def _delete_unknown_terms_from_ontology_graph(self, node: OntologyTree):
        """
        Deletes nodes that have no name in the CL ontology from the ontology graph.

        Arguments
        ---------
        node: dict
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        """
        new_children = []
        children = [] if node.children is None else node.children
        for child in children:
            unknown = child.name.startswith("CL:")
            if not unknown:
                new_children.append(child)
        if len(new_children) > 0:
            node.children = new_children
        else:
            node.children = None

        children = [] if node.children is None else node.children
        for child in children:
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
                subtree=OntologyTree(
                    id=node.id + "__" + str(node_count),
                    name=self.id_to_name[node.id] if node.id in self.id_to_name else node.id,
                    n_cells_rollup=int(
                        self.cell_counts_df_rollup[node.id] if node.id in self.cell_counts_df_rollup else 0
                    ),
                    n_cells=int(self.cell_counts_df[node.id] if node.id in self.cell_counts_df else 0),
                ),
                traverse_node_counter=None,
                all_unique_nodes=None,
            )

        children = []
        # sort subclasses such that if a node is in TRAVERSAL_PRIORITY_NODES, it is placed first
        subclasses.sort(key=lambda x: x.id not in TRAVERSAL_PRIORITY_NODES)

        for child in subclasses:
            traverse_ontology_result = self._traverse_ontology_with_counting(
                child, traverse_node_counter, all_unique_nodes
            )
            children.append(traverse_ontology_result.subtree)

        return TraverseOntologyResult(
            subtree=OntologyTree(
                id=node.id + "__" + str(node_count),
                name=self.id_to_name[node.id] if node.id in self.id_to_name else node.id,
                n_cells_rollup=int(self.cell_counts_df_rollup[node.id] if node.id in self.cell_counts_df_rollup else 0),
                n_cells=int(self.cell_counts_df[node.id] if node.id in self.cell_counts_df else 0),
                children=children,
            ),
            traverse_node_counter=traverse_node_counter,
            all_unique_nodes=all_unique_nodes,
        )

    def get_ontology_tree(self) -> OntologyTree:
        return self.ontology_graph

    ### Get the ontology tree state
    def get_ontology_tree_state_per_celltype(self) -> Dict[str, OntologyTreeState]:
        """
        This function gets the ontology tree state per cell type. The ontology tree state is a mask that determines
        which nodes are expanded by default and which nodes are not shown when expanded in the ontology.

        The following rules determine which cell types are expanded by default and not shown when expanded:
         - All instances of the target cell type are shown.
         - All children of the first instance of the target cell type are shown.
         - All grandchildren of the first instance of the target cell type are shown.
         - All siblings of the first instance of the target cell type are shown.
         - All nodes in the paths to any instance of the target cell type are shown.
        Nodes that are not shown when their parents are expanded will be collapsed into a dummy node that contains the
        ids of the children that are not shown. Nodes will be collapsed into dummy nodes if they are siblings of nodes
        that are in a path to any instances of the target cell type.

        Returns
        -------
        all_states_per_cell_type: dict
            A dictionary that maps cell type ontology term ids to their ontology tree states.
            The ontology tree state is a dictionary with keys
             - "isExpandedNodes": A list of cell type ontology term ids that are expanded by default.
             - "notShownWhenExpandedNodes": A dictionary that maps cell type ontology term ids to their hidden children
        """
        results = {
            celltype: delayed(self._process_cell_type__parallel)(celltype)
            for celltype in self.all_cell_type_ids_in_corpus
            if celltype in self.traverse_node_counter
        }

        logger.info(
            f"Getting ontology tree states for {len(results)} cell types using {CELLGUIDE_PIPELINE_NUM_CPUS} CPUs..."
        )
        with ProgressBar():
            computed_results = compute(*list(results.values()), num_workers=CELLGUIDE_PIPELINE_NUM_CPUS)
            all_states_per_cell_type = {
                cell_type_id: result for cell_type_id, result in zip(results.keys(), computed_results)
            }

        return all_states_per_cell_type

    def _process_cell_type__parallel(self, celltype: str) -> OntologyTreeState:
        """
        This function processes a cell type in parallel to generate its ontology tree state. The ontology tree state is a mask that determines
        which nodes are expanded by default and which nodes are not shown when expanded in the ontology.

        Arguments
        ---------
        celltype: str
            The cell type ontology term id.

        Returns
        -------
        OntologyTreeState
            The ontology tree state for the cell type. The ontology tree state is a dictionary with keys
                - "isExpandedNodes": A list of cell type ontology term ids that are expanded by default.
                - "notShownWhenExpandedNodes": A dictionary that maps cell type ontology term ids to their hidden children
        """

        all_paths = []
        siblings = []
        for i in range(self.traverse_node_counter[celltype]):
            end_node_i = celltype + "__" + str(i)
            path = self._depth_first_search_pathfinder(end_node_i)
            path = path if path else [end_node_i]
            all_paths.append(path)

            siblings.append(
                sum([self.all_children.get(parent, []) for parent in self.all_parents.get(end_node_i, [])], [])
            )

        visited_nodes_in_paths = list(set(sum(all_paths, [])))  # in a path to target

        children = self.all_children.get(celltype + "__0", [])  # children
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

        return OntologyTreeState(
            isExpandedNodes=isExpandedNodes,
            notShownWhenExpandedNodes=notShownWhenExpanded,
        )

    def get_ontology_tree_state_per_tissue(self) -> Dict[str, OntologyTreeState]:
        """
        This function gets the ontology tree state per tissue. The ontology tree state is a mask that determines
        which nodes are expanded by default and which nodes are not shown when expanded in the ontology.

        The following rules determine which cell types are expanded by default and not shown when expanded:
         - Only the first instance of each cell type present in the tissue is shown.
         - All nodes in a path to any of the first instances of the cell types present in the tissue are shown.
         - Nodes that are not outlier branches (depth=1 and the node has less than 0.1% observed cells out of the total number of cells) are shown.
        Nodes that are not shown when expanded will be collapsed into a dummy node that contains the ids of the
        children that are not shown. Nodes will be collapsed into dummy nodes if they are siblings of nodes that are
        in a path to the first instances of cell types present in the tissue.

        Returns
        -------
        all_states_per_tissue: dict
            A dictionary that maps tissue ontology term ids to their ontology tree states.
            The ontology tree state is a dictionary with keys
                - "isExpandedNodes": A list of cell type ontology term ids that are expanded by default.
                - "notShownWhenExpandedNodes": A dictionary that maps cell type ontology term ids to their hidden children
                - "tissueCounts": A dictionary that maps cell type ontology term ids to the number of cells in the tissue.
                    The number of cells is a dictionary with keys "n_cells" and "n_cells_rollup".
        """

        results = {
            tissue: delayed(self._process_tissue__parallel)(tissue)
            for tissue in self.uberon_by_celltype
            if tissue.startswith("UBERON:") and " (" not in tissue
        }

        logger.info(
            f"Getting ontology tree states for {len(results)} tissues using {CELLGUIDE_PIPELINE_NUM_CPUS} CPUs..."
        )
        with ProgressBar():
            computed_results = compute(*list(results.values()), num_workers=CELLGUIDE_PIPELINE_NUM_CPUS)
            all_states_per_tissue = {tissue_id: result for tissue_id, result in zip(results.keys(), computed_results)}

        return all_states_per_tissue

    def _process_tissue__parallel(self, tissueId: str) -> OntologyTreeState:
        """
        This function processes a tissue in parallel to generate its ontology tree state. The ontology tree state is a mask that determines
        which nodes are expanded by default and which nodes are not shown when expanded in the ontology.

        Arguments
        ---------
        tissueId: str
            The tissue ontology term id.

        Returns
        -------
        OntologyTreeState
            The ontology tree state for the tissue. The ontology tree state is a dictionary with keys
                - "isExpandedNodes": A list of cell type ontology term ids that are expanded by default.
                - "notShownWhenExpandedNodes": A dictionary that maps cell type ontology term ids to their hidden children
                - "tissueCounts": A dictionary that maps cell type ontology term ids to the number of cells in the tissue.
                    The number of cells is a dictionary with keys "n_cells" and "n_cells_rollup".
        """

        tissue_term = self.uberon_ontology[tissueId]
        tissue_label = tissue_term.name

        end_nodes = self.uberon_by_celltype[tissueId]
        uberon_ancestors = [i.id for i in tissue_term.superclasses()]

        # filter out hemaotoietic cell types from non-whitelisted tissues
        uberon_ancestors_in_whitelist = list(set(TISSUES_IMMUNE_CELL_WHITELIST).intersection(uberon_ancestors))
        if len(uberon_ancestors_in_whitelist) == 0:
            end_nodes_that_are_not_hematopoietic = [
                e
                for e in end_nodes
                if HEMATOPOIETIC_CELL_TYPE_ID not in [i.id for i in self.ontology[e].superclasses()]
            ]
            if len(end_nodes_that_are_not_hematopoietic) == 0:
                logger.info(f"Not filtering out immune cell for {tissue_label}")
            else:
                end_nodes = end_nodes_that_are_not_hematopoietic
        else:
            logger.info(f"Not filtering out immune cell for {tissue_label}")

        tissue_ct_df = self.tissue_celltypes_df[self.tissue_celltypes_df["tissue_ontology_term_id"] == tissueId]

        # attach cells not in the tissue cell types dataframe
        # these cells will be added with 0 counts
        df = tissue_ct_df[["cell_type_ontology_term_id", "n_cells"]]
        to_attach = pd.DataFrame()
        to_attach["cell_type_ontology_term_id"] = [
            i for i in self.all_cell_type_ids if i not in df["cell_type_ontology_term_id"].values
        ]
        to_attach["n_cells"] = 0
        df = pd.concat([df, to_attach], axis=0)
        # rollup the cell counts
        df["n_cells_rollup"] = df["n_cells"]
        df_rollup = rollup_across_cell_type_descendants(df, parallel=False, ignore_cols=["n_cells"])
        df_rollup = df_rollup[df_rollup["n_cells_rollup"] > 0]

        celltype_counts_in_tissue = dict(
            zip(
                df_rollup["cell_type_ontology_term_id"],
                df_rollup[["n_cells", "n_cells_rollup"]].to_dict(orient="records"),
            )
        )

        all_paths = []
        for end_node in end_nodes:
            if end_node in self.traverse_node_counter:
                end_node_0 = end_node + "__0"  # only get path to the first instance of a node.
                path = self._depth_first_search_pathfinder(end_node_0)
                path = path if path else [end_node_0]
                all_paths.append(path)

        valid_nodes = list(set(sum(all_paths, [])))

        ontology_graph_copy = self._get_deepcopy_of_ontology_graph()

        self._truncate_graph_in_tissue(
            graph=ontology_graph_copy,
            valid_nodes=valid_nodes,
            total_count=self.tissue_counts_df[tissueId],
            tissue_cell_counts=celltype_counts_in_tissue,
        )

        isExpandedNodes = list(set(_getExpandedData(ontology_graph_copy)))
        notShownWhenExpandedNodes = _getShownData(ontology_graph_copy)

        notShownWhenExpanded = {}
        for i in notShownWhenExpandedNodes:
            notShownWhenExpanded.update(i)

        return OntologyTreeState(
            isExpandedNodes=isExpandedNodes,
            notShownWhenExpandedNodes=notShownWhenExpanded,
            tissueCounts=celltype_counts_in_tissue,
        )

    def _depth_first_search_pathfinder(self, path_end_node, node=None, path=None) -> list[str]:
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

    def _truncate_graph_first_pass(self, graph: OntologyTree, valid_nodes: list[str]) -> bool:
        """
        This function truncates the ontology graph in-place by removing nodes that are not in the valid nodes list.
        It also adds a dummy node to the graph if the graph has children that are not in the valid nodes list.
        The dummy node contains the ids of the invalid children.

        Arguments
        ---------
        graph: OntologyTree
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        valid_nodes: list
            A list of valid cell type ontology term ids.

        Returns
        -------
        bool
            Returns True if the node is in valid nodes and False otherwise.
            This is only used for the recursion.
        """
        if graph.id not in valid_nodes:
            return False

        children = [] if graph.children is None else graph.children
        valid_children = []

        invalid_children_ids = []
        for child in children:
            is_valid = self._truncate_graph_first_pass(child, valid_nodes)
            if is_valid:
                valid_children.append(child)
            elif child.id != "":
                invalid_children_ids.append(child.id)

        if len(invalid_children_ids) > 0 and len(valid_children) > 0:
            valid_children.append(
                OntologyTree(
                    id="",
                    name="",
                    n_cells_rollup=0,
                    n_cells=0,
                    invalid_children_ids=invalid_children_ids,
                    parent=graph.id,
                )
            )
        if len(valid_children) > 0:
            graph.children = valid_children
        else:
            graph.children = None

        return True

    def _truncate_graph_second_pass(
        self, graph: OntologyTree, visited_nodes_in_paths: list[str], nodesWithChildrenFound: set[str] = None
    ) -> None:
        """
        We do not want to show cell types multiple times in the ontology tree unless they are nodes that are in one of the
        paths to the target cell type. In this case, show the node in the path and collapse all its siblings into the
        dummy node.

        Arguments
        ---------
        graph: OntologyTree
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
        if graph.id.split("__")[0] in nodesWithChildrenFound:
            # if node is in path to target cell type, preserve it
            # otherwise, collapse its siblings into a dummy node
            if graph.children is not None:
                children = graph.children
                new_children = []
                invalid_children_ids = []
                for child in children:
                    if child.id in visited_nodes_in_paths:
                        new_children.append(child)
                    elif child.id != "":
                        invalid_children_ids.append(child.id)

                if len(children) > len(new_children) and len(new_children) > 0:
                    # append dummy
                    new_children.append(
                        OntologyTree(
                            id="",
                            name="",
                            n_cells=0,
                            n_cells_rollup=0,
                            invalid_children_ids=invalid_children_ids,
                            parent=graph.id,
                        )
                    )
                if len(new_children) > 0:
                    graph.children = new_children
                else:
                    graph.children = None
        elif graph.children is not None:
            nodesWithChildrenFound.add(graph.id.split("__")[0])

        children = [] if graph.children is None else graph.children
        for child in children:
            if child.id != "":
                self._truncate_graph_second_pass(child, visited_nodes_in_paths, nodesWithChildrenFound)

    def _truncate_graph_in_tissue(
        self,
        *,
        graph: OntologyTree,
        valid_nodes: list[str],
        total_count: int,
        tissue_cell_counts: Dict[str, Dict[str, int]],
        seen_nodes_per_tissue: set[str] = None,
        depth: int = 0,
    ):
        """
        This function truncates the ontology graph in-place by removing nodes with any of the following properties:
         - The node is not in the valid nodes list.
         - The node is an outlier branch (depth=1 and the node has less than 0.1% observed cells out of the total number of cells).
         - The node is a child of a node that has already been encountered during the traversal.
        This function is used for the ontology tree state per tissue.

        Arguments
        ---------
        graph: OntologyTree
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        valid_nodes: list
            A list of valid cell type ontology term ids. These are cell types that are present in the tissue based on
            the CELLxGENE corpus.
        total_count: int
            The total number of cells in the tissue.
        tissue_cell_counts: dict
            A dictionary that maps cell type ontology term ids to the number of cells in the tissue.
        seen_nodes_per_tissue: set, optional, default=None
            A set of cell type ontology term ids that have been encountered during the traversal.
        depth: int, optional, default=0
            The depth of the current node in the ontology graph.
        """
        if seen_nodes_per_tissue is None:
            seen_nodes_per_tissue = set()

        children = [] if graph.children is None else graph.children
        if len(children):
            new_children = []
            invalid_children_ids = []
            for child in children:
                # if depth is 1 and the child node has less than 0.1% observed cells out of the total number of cells, then it is an outlier branch and must be pruned
                outlier_branch = (
                    depth == 1
                    and (
                        tissue_cell_counts.get(child.id.split("__")[0], {"n_cells_rollup": 0})["n_cells_rollup"]
                        / total_count
                        * 100
                    )
                    < 0.1
                )
                if child.id in valid_nodes and child.id not in seen_nodes_per_tissue and not outlier_branch:
                    new_children.append(child)
                    seen_nodes_per_tissue.add(child.id)
                else:
                    invalid_children_ids.append(child.id)
            if len(new_children) == 0:
                graph.children = None
            elif len(invalid_children_ids) > 0:
                new_children.append(
                    OntologyTree(
                        id="",
                        name="",
                        n_cells=0,
                        n_cells_rollup=0,
                        invalid_children_ids=invalid_children_ids,
                        parent=graph.id,
                    )
                )
                graph.children = new_children
            else:
                graph.children = new_children

            children = [] if graph.children is None else graph.children
            for child in children:
                if child.id != "":
                    self._truncate_graph_in_tissue(
                        graph=child,
                        valid_nodes=valid_nodes,
                        total_count=total_count,
                        tissue_cell_counts=tissue_cell_counts,
                        seen_nodes_per_tissue=seen_nodes_per_tissue,
                        depth=depth + 1,
                    )

    def _initialize_children_and_parents_per_node(self) -> Tuple[Dict[str, list[str]], Dict[str, list[str]]]:
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

    def _build_children_per_node(self, node: OntologyTree, all_children=None) -> Dict[str, list[str]]:
        """
        Builds a dictionary that maps cell type ontology term ids to their children.

        Arguments
        ---------
        node: OntologyTree
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
        children = [] if node.children is None else node.children
        ids = [] if len(children) == 0 else [child.id for child in children]

        all_children[node.id] = ids
        for child in children:
            self._build_children_per_node(child, all_children)
        return all_children

    def _build_parents_per_node(self, node: OntologyTree, all_parents=None) -> Dict[str, list[str]]:
        """
        Builds a dictionary that maps cell type ontology term ids to their parents.

        Arguments
        ---------
        node: OntologyTree
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
        children = [] if node.children is None else node.children

        for child in children:
            all_parents[child.id] = [node.id]
            self._build_parents_per_node(child, all_parents)
        return all_parents

    def _get_deepcopy_of_ontology_graph(self) -> OntologyTree:
        """
        Returns a deepcopy of the ontology graph.
        """
        return copy.deepcopy(self.ontology_graph)


def _getExpandedData(ontology_graph: OntologyTree, isExpandedNodes=None) -> list[str]:
    """
    This function gets the cell type ontology term ids of the nodes that are expanded by default.

    Arguments
    ---------
    ontology_graph: OntologyTree
        A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.

    Returns
    -------
    isExpandedNodes: list
        A list of cell type ontology term ids that are expanded by default.
    """
    if isExpandedNodes is None:
        isExpandedNodes = []

    if ontology_graph.children is not None:
        isExpandedNodes.append(ontology_graph.id)
        for child in ontology_graph.children:
            _getExpandedData(child, isExpandedNodes)

    return isExpandedNodes


def _getShownData(graph: OntologyTree, notShownWhenExpandedNodes=None) -> Dict[str, Dict[str, list[str]]]:
    """
    This function gets the children nodes that are not shown when expanded for each node in the ontology graph.

    Arguments
    ---------
    graph: OntologyTree
        A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.

    Returns
    -------
    notShownWhenExpandedNodes: list
        A list of dictionaries that map cell type ontology term ids to their hidden children.
    """
    if notShownWhenExpandedNodes is None:
        notShownWhenExpandedNodes = []

    if graph.children is not None:
        for child in graph.children:
            if child.id == "":
                if child.invalid_children_ids is not None and len(child.invalid_children_ids) > 0:
                    notShownWhenExpandedNodes.append({child.parent: list(set(child.invalid_children_ids))})
            else:
                _getShownData(child, notShownWhenExpandedNodes)
    return notShownWhenExpandedNodes


def _to_dict(a: list[Any], b: list[Any]) -> Dict[Any, list[Any]]:
    """
    convert a flat key array (a) and a value array (b) into a dictionary with values grouped by keys
    """
    a = np.array(a)
    b = np.array(b)
    idx = np.argsort(a)
    a = a[idx]
    b = b[idx]
    bounds = np.where(a[:-1] != a[1:])[0] + 1
    bounds = np.append(np.append(0, bounds), a.size)
    bounds_left = bounds[:-1]
    bounds_right = bounds[1:]
    slists = [b[bounds_left[i] : bounds_right[i]] for i in range(bounds_left.size)]
    d = dict(zip(np.unique(a), [list(set(x)) for x in slists]))
    return d
