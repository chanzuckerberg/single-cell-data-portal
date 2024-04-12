from dataclasses import dataclass
from typing import Optional


@dataclass
class GoTermMetadata:
    name: str
    id: str
    goDescription: Optional[str]
    synonyms: list[str]
