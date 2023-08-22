from dataclasses import dataclass, field
from typing import Dict, List, Optional

from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder


@dataclass
class OntologyTreeState:
    isExpandedNodes: list[str]
    notShownWhenExpandedNodes: Dict[str, Dict[str, list[str]]]
    tissueCounts: Optional[Dict[str, int]] = field(default=None)


@dataclass
class OntologyTree:
    id: str
    name: str
    n_cells_rollup: int
    n_cells: int
    children: Optional[List["OntologyTree"]] = field(default=None)
    invalid_children_ids: Optional[List[str]] = field(default=None)
    parent: Optional[str] = field(default=None)


@dataclass
class OntologyTreeData:
    tree_builder: OntologyTreeBuilder
    ontology_graph: OntologyTree
    all_states_per_cell_type: OntologyTreeState
    all_states_per_tissue: OntologyTreeState
