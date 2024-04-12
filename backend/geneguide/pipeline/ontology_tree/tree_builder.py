import copy
import logging
from collections import deque
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from dask import compute, delayed
from dask.diagnostics import ProgressBar
from pronto import Ontology, Term

from backend.geneguide.pipeline.constants import GENEGUIDE_PIPELINE_NUM_CPUS, GO_URL
from backend.geneguide.pipeline.ontology_tree.types import OntologyTree, OntologyTreeState

logger = logging.getLogger(__name__)


@dataclass
class TraverseOntologyResult:
    subtree: OntologyTree
    traverse_node_counter: Dict[str, int]
    all_unique_nodes: set[str]


GO_ROOT_NODES = [
    "GO:0008150",  # biological_process
    "GO:0005575",  # cellular_component
    "GO:0003674",  # molecular_function
]


class OntologyTreeBuilder:
    def __init__(self):

        logger.info("Loading GO ontology...")
        self.ontology = Ontology(GO_URL)
        self.id_to_name = self._initialize_id_to_name()
        self.all_go_ids = list(self.id_to_name.keys())

        logger.info("Initializing ontology tree data structures by traversing GO ontology...")

        self.all_unique_nodes = set()
        self.traverse_node_counter = {}
        children = []
        for root_node in GO_ROOT_NODES:
            traverse_ontology_result = self._traverse_ontology_with_counting(
                self.ontology[root_node],
                traverse_node_counter=self.traverse_node_counter,
                all_unique_nodes=self.all_unique_nodes,
            )
            children.append(traverse_ontology_result.subtree)

        self.ontology_graph = OntologyTree(id="GO:0000000__0", name="Root", children=children)

        self.all_children, self.all_parents = self._initialize_children_and_parents_per_node()
        self.start_node = "GO:0000000__0"

    def _initialize_id_to_name(self):
        """
        Initializes a dictionary that maps GO term ontology term ids to their names.

        Returns
        -------
        id_to_name: dict
            A dictionary that maps GO term ontology term ids to their names.
        all_cell_type_ids: list
            A list of all GO term ontology term ids.
        """
        id_to_name = {}
        for term in self.ontology.terms():
            if term.id.startswith("GO:"):
                id_to_name[term.id] = term.name
        return id_to_name

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
            A dictionary that maps GO term ontology term ids to the number of times they have been visited.
        all_unique_nodes: set, optional, default=None
            A set of all unique GO term ontology term ids that have been visited.

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
                ),
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
            subtree=OntologyTree(
                id=node.id + "__" + str(node_count),
                name=self.id_to_name[node.id] if node.id in self.id_to_name else node.id,
                children=children,
            ),
            traverse_node_counter=traverse_node_counter,
            all_unique_nodes=all_unique_nodes,
        )

    def get_ontology_tree(self) -> OntologyTree:
        return self.ontology_graph

    ### Get the ontology tree state
    def get_ontology_tree_state_per_goterm(self) -> Dict[str, OntologyTreeState]:
        """
        This function gets the ontology tree state per GO term. The ontology tree state is a mask that determines
        which nodes are expanded by default and which nodes are not shown when expanded in the ontology.

        The following rules determine which GO terms are expanded by default and not shown when expanded:
         - All instances of the target GO term are shown.
         - All children of the first instance of the target GO term are shown.
         - All grandchildren of the first instance of the target GO term are shown.
         - All siblings of the first instance of the target GO term are shown.
         - All nodes in the paths to any instance of the target GO term are shown.
        Nodes that are not shown when their parents are expanded will be collapsed into a dummy node that contains the
        ids of the children that are not shown. Nodes will be collapsed into dummy nodes if they are siblings of nodes
        that are in a path to any instances of the target GO term.

        Returns
        -------
        all_states_per_go_term: dict
            A dictionary that maps GO term ontology term ids to their ontology tree states.
            The ontology tree state is a dictionary with keys
             - "isExpandedNodes": A list of GO term ontology term ids that are expanded by default.
             - "notShownWhenExpandedNodes": A dictionary that maps GO term ontology term ids to their hidden children
        """
        results = {goterm: delayed(self._process_go_term__parallel)(goterm) for goterm in self.traverse_node_counter}

        logger.info(
            f"Getting ontology tree states for {len(results)} GO terms using {GENEGUIDE_PIPELINE_NUM_CPUS} CPUs..."
        )
        with ProgressBar():
            computed_results = compute(*list(results.values()), num_workers=GENEGUIDE_PIPELINE_NUM_CPUS)
            all_states_per_goterm = {goterm_id: result for goterm_id, result in zip(results.keys(), computed_results)}

        return all_states_per_goterm

    def _process_go_term__parallel(self, goterm: str) -> OntologyTreeState:
        """
        This function processes a GO term in parallel to generate its ontology tree state. The ontology tree state is a mask that determines
        which nodes are expanded by default and which nodes are not shown when expanded in the ontology.

        Arguments
        ---------
        goterm: str
            The GO term id.

        Returns
        -------
        OntologyTreeState
            The ontology tree state for the GO term. The ontology tree state is a dictionary with keys
                - "isExpandedNodes": A list of GO term ontology term ids that are expanded by default.
                - "notShownWhenExpandedNodes": A dictionary that maps GO term ontology term ids to their hidden children
        """

        all_paths = []
        siblings = []
        for i in range(self.traverse_node_counter[goterm]):
            end_node_i = goterm + "__" + str(i)
            path = self._depth_first_search_pathfinder(end_node_i)
            path = path if path else [end_node_i]
            all_paths.append(path)

            siblings.append(
                sum([self.all_children.get(parent, []) for parent in self.all_parents.get(end_node_i, [])], [])
            )

        visited_nodes_in_paths = list(set(sum(all_paths, [])))  # in a path to target

        children = self.all_children.get(goterm + "__0", [])  # children
        grandchildren = sum([self.all_children.get(child, []) for child in children], [])  # grandchildren
        siblings = list(set(sum(siblings, [])))  # siblings
        valid_nodes = list(set(visited_nodes_in_paths + children + grandchildren + siblings))

        new_ontology_graph = self._build_truncated_graph_first_pass(valid_nodes)
        self._truncate_graph_second_pass(new_ontology_graph, visited_nodes_in_paths)

        isExpandedNodes = list(set(_getExpandedData(new_ontology_graph)))
        notShownWhenExpandedNodes = _getShownData(new_ontology_graph)

        notShownWhenExpanded = {}
        for i in notShownWhenExpandedNodes:
            notShownWhenExpanded.update(i)

        return OntologyTreeState(
            isExpandedNodes=isExpandedNodes,
            notShownWhenExpandedNodes=notShownWhenExpanded,
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

    def _build_truncated_graph_first_pass(self, valid_nodes: set[str]) -> Optional[OntologyTree]:

        # instantiate a new graph
        new_graph = OntologyTree(id=self.ontology_graph.id, name=self.ontology_graph.name, children=[])
        queue = deque([(self.ontology_graph, new_graph)])

        while queue:
            # pop the node off the queue
            node, new_node = queue.popleft()

            # if the node is not in the valid nodes list, skip it
            if node.id not in valid_nodes:
                continue

            # get the children of the node
            children = [] if node.children is None else node.children
            invalid_children_ids = []

            # for each child
            for child in children:
                # create a new child node
                new_child = OntologyTree(id=child.id, name=child.name, children=[])

                # if the child is in the valid nodes list, add it to the queue
                if child.id in valid_nodes:
                    queue.append((child, new_child))

                    # add the child to the new node
                    new_node.children.append(new_child)
                elif child.id != "":
                    # if the child is not in the valid nodes list, add its id to the invalid children ids
                    invalid_children_ids.append(child.id)

            # if the new node has invalid children, add a dummy node to the new node
            if len(invalid_children_ids) > 0 and len(new_node.children) > 0:
                new_node.children.append(
                    OntologyTree(
                        id="",
                        name="",
                        invalid_children_ids=invalid_children_ids,
                        parent=node.id,
                    )
                )
        return new_graph

    def _truncate_graph_second_pass(
        self, graph: OntologyTree, visited_nodes_in_paths: list[str], nodesWithChildrenFound: set[str] = None
    ) -> None:
        """
        We do not want to show GO terms multiple times in the ontology tree unless they are nodes that are in one of the
        paths to the target GO term. In this case, show the node in the path and collapse all its siblings into the
        dummy node.

        Arguments
        ---------
        graph: OntologyTree
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        visited_nodes_in_paths: list
            A list of GO term ontology term ids that are in one of the paths to the target GO term.
            The target GO term may have multiple instances in the ontology graph (each instance with a unique suffix).
        nodesWithChildrenFound: set, optional, default=None
            A set of GO term ontology term ids that have children. This is used to remove children of nodes that
            have already been encountered during the traversal.
        """
        if nodesWithChildrenFound is None:
            nodesWithChildrenFound = set()
        if graph.id.split("__")[0] in nodesWithChildrenFound:
            # if node is in path to target GO term, preserve it
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

    def _build_truncated_graph_second_pass(
        self, node: OntologyTree, visited_nodes_in_paths: set[str], parent_id: str = None
    ) -> Optional[OntologyTree]:
        if node.id not in visited_nodes_in_paths:
            return None

        new_node = OntologyTree(id=node.id, name=node.name, children=[])

        for child in node.children:
            refined_child = self._refine_graph_second_pass(child, visited_nodes_in_paths, parent_id=node.id)
            if refined_child:
                new_node.children.append(refined_child)
            elif child.id and child.id not in visited_nodes_in_paths:
                # Handle dummy node creation or other specific logic here
                pass

        return new_node

    def _initialize_children_and_parents_per_node(self) -> Tuple[Dict[str, list[str]], Dict[str, list[str]]]:
        """
        Initializes dictionaries that map GO term ontology term ids to their children and parents.
        While the gene ontology contains multiple instances of the same node, the ontology graph
        contains only one instance of each node as each GO term has been suffixed with a unique integer.
        Therefore, each node only has one parent and a unique set of children.

        Returns
        -------
        all_children: dict
            A dictionary that maps GO term ontology term ids to their children.
        all_parents: dict
            A dictionary that maps GO term ontology term ids to their parents.
        """
        all_children = self._build_children_per_node(self.ontology_graph)
        all_parents = self._build_parents_per_node(self.ontology_graph)
        return all_children, all_parents

    def _build_children_per_node(self, node: OntologyTree, all_children=None) -> Dict[str, list[str]]:
        """
        Builds a dictionary that maps GO term ontology term ids to their children.

        Arguments
        ---------
        node: OntologyTree
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        all_children: dict, optional, default=None
            A dictionary that maps GO term ontology term ids to their children.

        Returns
        -------
        all_children: dict
            A dictionary that maps GO term ontology term ids to their children.
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
        Builds a dictionary that maps GO term ontology term ids to their parents.

        Arguments
        ---------
        node: OntologyTree
            A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.
        all_parents: dict, optional, default=None
            A dictionary that maps GO term ontology term ids to their parents.

        Returns
        -------
        all_parents: dict
            A dictionary that maps GO term ontology term ids to their parents.
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
    This function gets the GO term ontology term ids of the nodes that are expanded by default.

    Arguments
    ---------
    ontology_graph: OntologyTree
        A node of the ontology graph. At the top-level, this should be the root node of the ontology graph.

    Returns
    -------
    isExpandedNodes: list
        A list of GO term ontology term ids that are expanded by default.
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
        A list of dictionaries that map GO term ontology term ids to their hidden children.
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
