from dataclasses import dataclass
from typing import Optional


@dataclass
class CellMetadata:
    name: str
    id: str
    clDescription: Optional[str]
    synonyms: list[str]
