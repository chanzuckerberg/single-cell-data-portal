from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class OntologyTreeStatePerCellType:
    isExpandedNodes: Dict[str, list[str]]
    notShownWhenExpandedNodes: Dict[str, Dict[str, list[str]]]


@dataclass
class OntologyTreeStatePerTissue:
    isExpandedNodes: Dict[str, list[str]]
    notShownWhenExpandedNodes: Dict[str, Dict[str, list[str]]]
    tissueCounts: Dict[str, int]


@dataclass
class OntologyTree:
    id: str
    name: str
    n_cells_rollup: int
    n_cells: int
    children: Optional[List["OntologyTree"]]
