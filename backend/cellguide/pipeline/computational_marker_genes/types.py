from dataclasses import dataclass
from typing import Dict


@dataclass
class ComputationalMarkerGenes:
    me: float
    pc: float
    marker_score: float
    symbol: str
    name: str
    groupby_dims: Dict[str, str]
