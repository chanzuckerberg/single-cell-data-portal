from dataclasses import dataclass, field
from typing import Optional

"""
This module defines the data classes used in the canonical marker genes pipeline.

All the data classes except for ParsedAsctbTableEntry correspond to the structure of the ASCTB data maintained
by HUBMAP. This is explicitly defined to raise errors should ASCTB update the structure of their master ASCTB JSON.

Classes:
    AnatomicalStructure: Represents an anatomical structure with an ID and label as structured in ASCTB.
    GeneBiomarker: Represents a gene biomarker with an ID, label, name, type, and optional notes as structured in ASCTB.
    Reference: Represents a reference with an ID, DOI, and notes as structured in ASCTB.
    ParsedAsctbTableEntry: Represents a parsed ASCTB table entry with tissue, symbol, name, publication, publication titles, and cell type ontology term ID.
"""


@dataclass
class AnatomicalStructure:
    id: str
    rdfs_label: str
    name: str
    notes: Optional[str] = field(default=None)


@dataclass
class GeneBiomarker:
    id: str
    rdfs_label: str
    name: str
    b_type: str
    notes: Optional[str] = field(default=None)


@dataclass
class Reference:
    id: str
    doi: str
    notes: Optional[str] = field(default=None)


@dataclass
class ParsedAsctbTableEntry:
    tissue: str
    symbol: str
    name: str
    publication: str
    publication_titles: str
