import json

import numpy as np
import pandas as pd
from pronto import Ontology

from backend.wmg.data.constants import CL_BASIC_PERMANENT_URL_PRONTO, UBERON_BASIC_PERMANENT_URL_PRONTO
from backend.wmg.data.rollup import rollup_across_cell_type_descendants

# tissues in which we want to show immune cell types
TISSUES_IMMUNE_CELL_WHITELIST = [
    "UBERON:0002193",  # hemolymphoid system
    "UBERON:0002390",  # hematopoietic system
    "UBERON:0000178",  # blood
    "UBERON:0005057",  # immune organ
]
# immune cell type ancestor to be filtered out from non-whitelisted tissues
HEMATOPOIETIC_CELL_TYPE_ID = "CL:0000988"


class OntologyTreeClimber:
    def __init__(self, cell_counts_df, root_node="CL:0000548"):
        self.ontology = Ontology(CL_BASIC_PERMANENT_URL_PRONTO)

        self.tissue_counts_df = cell_counts_df.groupby("tissue_ontology_term_id").sum(numeric_only=True)["n_cells"]
        self.tissue_celltypes_df = (
            cell_counts_df.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"])
            .sum(numeric_only=True)
            .reset_index()
        )

        self.uberon_by_celltype = _to_dict(
            cell_counts_df["tissue_ontology_term_id"], cell_counts_df["cell_type_ontology_term_id"]
        )

        self.id_to_name, self.all_cell_type_ids = self._initialze_id_to_name()
        self.cell_counts_df, self.cell_counts_df_rollup = self._initialize_cell_counts_df_rollup(cell_counts_df)

        self.ontology_graph, self.traverse_node_counter, self.all_unique_nodes = self._traverse_ontology_with_counting(
            self.ontology[root_node]
        )

        self._delete_unknown_terms_from_ontology_graph(self.ontology_graph)

        self.all_children, self.all_parents = self._initialize_children_and_parents_per_node()

        self.start_node = f"{root_node}__0"

    def _initialize_uberon_ontology(self):
        self.uberon_ontology = Ontology(UBERON_BASIC_PERMANENT_URL_PRONTO)
        self._initialze_id_to_name()

    def _initialize_cell_counts_df_rollup(self, cell_counts_df):
        cell_counts_df = cell_counts_df.groupby("cell_type_ontology_term_id").sum(numeric_only=True)[["n_cells"]]

        to_attach = pd.DataFrame()
        to_attach["cell_type_ontology_term_id"] = [
            i for i in self.all_cell_types_ids if i not in cell_counts_df["cell_type_ontology_term_id"]
        ]
        to_attach["n_cells"] = 0
        cell_counts_df = pd.concat([cell_counts_df, to_attach], axis=0)
        cell_counts_df_rollup = rollup_across_cell_type_descendants(self.cell_counts_df)
        cell_counts_df = cell_counts_df.set_index("cell_type_ontology_term_id")["n_cells"]
        cell_counts_df_rollup = cell_counts_df_rollup.set_index("cell_type_ontology_term_id")["n_cells"]
        return cell_counts_df, cell_counts_df_rollup

    def _get_deepcopy_of_ontology_graph(self):
        return json.loads(json.dumps(self.ontology_graph))

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
        return id_to_name, list(self.id_to_name)

    def _initialize_children_and_parents_per_node(self):
        all_children = self._build_children_per_node(self.ontology_graph)
        all_parents = self._build_parents_per_node(self.ontology_graph)
        return all_children, all_parents

    def _build_children_per_node(self, node, all_children=None):
        if all_children is None:
            all_children = {}
        children = node.get("children", [])
        ids = [] if len(children) == 0 else [child["id"] for child in children]

        all_children[node["id"]] = ids
        for child in children:
            self._build_children_per_node(child, all_children)
        return all_children

    def _build_parents_per_node(self, node, all_parents=None):
        if all_parents is None:
            all_parents = {}
        children = node.get("children", [])

        for child in children:
            all_parents[child["id"]] = [node["id"]]
            self._build_parents_per_node(child, all_parents)
        return all_parents

    def _traverse_ontology_with_counting(self, node, traverse_node_counter=None, all_unique_nodes=None):
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

    def _truncate_graph_second_pass(self, graph, visited_nodes_in_paths, nodesWithChildrenFound=None):
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

    def _truncate_graph_in_tissue(
        self, graph, valid_nodes, total_count, tissue_cell_counts, seen_nodes_per_tissue=None, depth=0
    ):
        if seen_nodes_per_tissue is None:
            seen_nodes_per_tissue = set()

        children = graph.get("children", [])
        if len(children):
            new_children = []
            invalid_children_ids = []
            for child in children:
                # if depth is 1 and the child node has less than 0.1% observed cells out of the total number of cells, then it is an outlier branch and must be pruned
                outlier_branch = (
                    depth == 1
                    and (
                        tissue_cell_counts.get(child["id"].split("__")[0], {"n_cells_rollup": 0})["n_cells_rollup"]
                        / total_count
                        * 100
                    )
                    < 0.1
                )
                if child["id"] in valid_nodes and child["id"] not in seen_nodes_per_tissue and not outlier_branch:
                    new_children.append(child)
                    seen_nodes_per_tissue.add(child["id"])
                else:
                    invalid_children_ids.append(child["id"])
            if len(new_children) == 0:
                del graph["children"]
            elif len(invalid_children_ids) > 0:
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
                graph["children"] = new_children
            else:
                graph["children"] = new_children

            for child in graph.get("children", []):
                if child["id"] != "":
                    self._truncate_graph_in_tissue(
                        child,
                        valid_nodes,
                        total_count,
                        tissue_cell_counts,
                        seen_nodes_per_tissue=seen_nodes_per_tissue,
                        depth=depth + 1,
                    )

    def get_ontology_tree_state_per_celltype(self):
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

    def get_ontology_tree_state_per_tissue(self):
        all_states_per_tissue = {}
        tissue_by_cell_type = []
        for tissue in self.uberon_by_celltype:
            if " (" not in tissue:
                print(tissue)
                tissueId = tissue
                tissue_term = self.uberon_ontology[tissueId]
                tissue_label = tissue_term.name

                end_nodes = self.uberon_by_celltype[tissue]
                uberon_ancestors = [i.id for i in tissue_term.superclasses()]
                if len(list(set(TISSUES_IMMUNE_CELL_WHITELIST).intersection(uberon_ancestors))) == 0:
                    end_nodes2 = [
                        e
                        for e in end_nodes
                        if HEMATOPOIETIC_CELL_TYPE_ID not in [i.id for i in self.ontology[e].superclasses()]
                    ]
                    if len(end_nodes2) == 0:
                        print("Not filtering out immune cell for", tissue_label)
                    else:
                        end_nodes = end_nodes2
                else:
                    print("Not filtering out immune cell for", tissue_label)

                tissue_ct_df = self.tissue_celltypes_df[self.tissue_celltypes_df["tissue_ontology_term_id"] == tissue]

                df = tissue_ct_df[["cell_type_ontology_term_id", "n_cells"]]
                to_attach = pd.DataFrame()
                to_attach["cell_type_ontology_term_id"] = [
                    i for i in self.all_cell_types_ids if i not in df["cell_type_ontology_term_id"].values
                ]
                to_attach["n_cells"] = 0
                df = pd.concat([df, to_attach], axis=0)
                df["n_cells_rollup"] = df["n_cells"]
                df_rollup = rollup_across_cell_type_descendants(df, ignore_cols=["n_cells"])
                df_rollup = df_rollup[df_rollup["n_cells_rollup"] > 0]

                celltype_counts_in_tissue = dict(
                    zip(
                        df_rollup["cell_type_ontology_term_id"],
                        df_rollup[["n_cells", "n_cells_rollup"]].to_dict(orient="records"),
                    )
                )

                tissue_by_cell_type.append({"id": tissue, "label": tissue_label})

                all_paths = []
                for end_node in end_nodes:
                    if end_node in self.traverse_node_counter:
                        end_node_0 = end_node + "__0"  # only get path to the first instance of a node.
                        path = self._depth_first_search_pathfinder(end_node_0)
                        path = [i[::-1] for i in path] if path else [end_node_0]
                        all_paths.append(path)

                valid_nodes = list(set(sum(all_paths, [])))

                ontology_graph_copy = self._get_deepcopy_of_ontology_graph()

                self._truncate_graph_in_tissue(
                    ontology_graph_copy, valid_nodes, self.tissue_counts_df[tissue], celltype_counts_in_tissue
                )

                isExpandedNodes = list(set(_getExpandedData(ontology_graph_copy)))
                notShownWhenExpandedNodes = _getShownData(ontology_graph_copy)

                notShownWhenExpanded = {}
                for i in notShownWhenExpandedNodes:
                    notShownWhenExpanded.update(i)

                all_states_per_tissue[tissue] = {
                    "isExpandedNodes": isExpandedNodes,
                    "notShownWhenExpandedNodes": notShownWhenExpanded,
                    "tissueCounts": celltype_counts_in_tissue,
                }
        return all_states_per_tissue


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


def _to_dict(a, b):
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
