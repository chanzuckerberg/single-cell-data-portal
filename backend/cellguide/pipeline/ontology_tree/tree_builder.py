import json

import pandas as pd
from pronto import Ontology

from backend.wmg.data.constants import CL_BASIC_PERMANENT_URL_PRONTO
from backend.wmg.data.rollup import rollup_across_cell_type_descendants


class OntologyTreeClimber:
    def __init__(self, cell_counts_df, root_node="CL:0000548"):
        self.ontology = Ontology(CL_BASIC_PERMANENT_URL_PRONTO)
        self._initialze_id_to_name()

        to_attach = pd.DataFrame()
        to_attach["cell_type_ontology_term_id"] = [
            i for i in self.all_cell_types_ids if i not in cell_counts_df["cell_type_ontology_term_id"]
        ]
        to_attach["n_cells"] = 0
        cell_counts_df = pd.concat([cell_counts_df, to_attach], axis=0)
        cell_counts_df_rollup = rollup_across_cell_type_descendants(self.cell_counts_df)
        cell_counts_df = cell_counts_df.set_index("cell_type_ontology_term_id")["n_cells"]
        cell_counts_df_rollup = cell_counts_df_rollup.set_index("cell_type_ontology_term_id")["n_cells"]

        self.cell_counts_df = cell_counts_df
        self.cell_counts_df_rollup = cell_counts_df_rollup

        self.traverse_node_counter = {}
        self.all_unique_nodes = set()
        self.ontology_graph = self._traverse_ontology_with_counting(self.ontology[root_node])
        self._delete_unknown_terms_from_ontology_graph(self.ontology_graph)

        self.all_children = {}
        self.all_parents = {}
        self._initialize_children_and_parents_per_node()

        self.start_node = f"{root_node}__0"

    def _delete_unknown_terms_from_ontology_graph(self, node):
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

    def _initialze_id_to_name(self):
        id_to_name = {}
        for term in self.ontology.terms():
            id_to_name[term.id] = term.name
        self.id_to_name = id_to_name
        self.all_cell_types_ids = list(self.id_to_name)

    def _initialize_children_and_parents_per_node(self):
        self._build_children_per_node(self.ontology_graph)
        self._build_parents_per_node(self.ontology_graph)

    def _build_children_per_node(self, node):
        children = node.get("children", [])
        ids = [] if len(children) == 0 else [child["id"] for child in children]

        self.all_children[node["id"]] = ids
        for child in children:
            self._build_children_per_node(child)

    def _build_parents_per_node(self, node):
        children = node.get("children", [])

        for child in children:
            self.all_parents[child["id"]] = [node["id"]]
            self._build_parents_per_node(child)

    def _traverse_ontology_with_counting(self, node):
        node_count = self.traverse_node_counter.get(node.id, 0)
        self.traverse_node_counter[node.id] = node_count + 1
        self.all_unique_nodes.add(node.id + "__" + str(node_count))

        subclasses = list(node.subclasses(with_self=False, distance=1))

        if len(subclasses) == 0:
            return {
                "id": node.id + "__" + str(node_count),
                "name": self.id_to_name[node.id] if node.id in self.id_to_name else node.id,
                "n_cells_rollup": int(
                    self.cell_counts_df_rollup[node.id] if node.id in self.cell_counts_df_rollup else 0
                ),
                "n_cells": int(self.cell_counts_df[node.id] if node.id in self.cell_counts_df else 0),
            }

        children = []
        for child in subclasses:
            children.append(self.traverse_ontology_with_counting(child))

        return {
            "id": node.id + "__" + str(node_count),
            "name": self.id_to_name[node.id] if node.id in self.id_to_name else node.id,
            "n_cells_rollup": int(self.cell_counts_df_rollup[node.id] if node.id in self.cell_counts_df_rollup else 0),
            "n_cells": int(self.cell_counts_df[node.id] if node.id in self.cell_counts_df else 0),
            "children": children,
        }

    def _depth_first_search_pathfinder(self, path_end_node, node=None, path=None):
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

    def _truncate_graph_second_pass(self, graph, visited_nodes_in_paths, nodesWithChildrenFound):
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

    def _get_deepcopy_of_ontology_graph(self):
        return json.loads(json.dumps(self.ontology_graph))

    def get_ontology_tree_state_per_celltype(self):
        # We only want to show terms that are CHILDREN, GRANDCHILDREN, SIBLINGS OF TARGET, or IN A PATH TO TARGET

        all_states_per_cell_type = {}
        for i, end_node in enumerate(self.all_cell_types_ids):
            if end_node in self.traverse_node_counter:
                all_paths = []
                siblings = []
                for i in range(self.traverse_node_counter[end_node]):
                    end_node_i = end_node + "__" + str(i)
                    path = self._depth_first_search_pathfinder(end_node_i)
                    path = [i[::-1] for i in path] if path else [end_node_i]
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
                self._truncate_graph_second_pass(ontology_graph_copy, visited_nodes_in_paths, set())
                self._delete_unknown_terms_from_ontology_graph(ontology_graph_copy)

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


def _getExpandedData(ontology_graph, isExpandedNodes=None):
    if isExpandedNodes is None:
        isExpandedNodes = []

    if "children" in ontology_graph:
        isExpandedNodes.append(ontology_graph["id"])
        for child in ontology_graph["children"]:
            _getExpandedData(child, isExpandedNodes)

    return isExpandedNodes


def _getShownData(graph, notShownWhenExpandedNodes=None):
    if notShownWhenExpandedNodes is None:
        notShownWhenExpandedNodes = []

    if "children" in graph:
        for child in graph["children"]:
            if child["id"] == "":
                if len(child["invalid_children_ids"]) > 0:
                    notShownWhenExpandedNodes.append({child["parent"]: list(set(child["invalid_children_ids"]))})
            else:
                _getShownData(child, notShownWhenExpandedNodes)
