from dataclasses import dataclass


@dataclass
class SourceCollectionsData:
    collection_name: str
    collection_url: str
    publication_url: str
    publication_title: str
    tissue: list[dict[str, str]]
    disease: list[dict[str, str]]
    organism: list[dict[str, str]]
