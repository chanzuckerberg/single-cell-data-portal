from dataclasses import dataclass
from typing import Optional


@dataclass
class CellMetadata:
    name: str
    id: str
    clDescription: Optional[str]
    synonyms: list[str]


@dataclass
class TissueMetadata:
    name: str
    id: str
    uberonDescription: Optional[str]
    synonyms: list[str]
