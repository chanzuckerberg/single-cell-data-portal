from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class OntologyTreeState:
    isExpandedNodes: Dict[str, list[str]]
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
