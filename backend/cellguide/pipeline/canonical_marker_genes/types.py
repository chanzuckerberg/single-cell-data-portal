from dataclasses import dataclass


@dataclass
class AnatomicalStructure:
    id: str
    rdfs_label: str
    name: str


@dataclass
class GeneBiomarker:
    id: str
    rdfs_label: str
    name: str
    b_type: str


@dataclass
class Reference:
    id: str
    doi: str
    notes: str


@dataclass
class ParsedAsctbTableEntry:
    tissue: str
    symbol: str
    name: str
    publication: str
    publication_titles: str
    cell_type_ontology_term_id: str
