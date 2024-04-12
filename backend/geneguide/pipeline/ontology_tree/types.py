from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class OntologyTreeState:
    isExpandedNodes: list[str]
    notShownWhenExpandedNodes: Dict[str, Dict[str, list[str]]]


@dataclass
class OntologyTree:
    id: str
    name: str
    children: Optional[List["OntologyTree"]] = field(default=None)
    invalid_children_ids: Optional[List[str]] = field(default=None)
    parent: Optional[str] = field(default=None)


@dataclass
class OntologyTreeData:
    ontology_graph: OntologyTree
    all_states_per_go_term: OntologyTreeState
