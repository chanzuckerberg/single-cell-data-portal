from dataclasses import dataclass


@dataclass
class ComputationalMarkerGenes:
    me: float
    pc: float
    marker_score: float
    symbol: str
    name: str
    groupby_dims: dict[str, str]
