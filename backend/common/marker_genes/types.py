from dataclasses import dataclass


@dataclass
class ComputationalMarkerGenes:
    me: float
    pc: float
    marker_score: float
    specificity: float
    gene_ontology_term_id: str
    symbol: str
    name: str
    groupby_dims: dict[str, str]
