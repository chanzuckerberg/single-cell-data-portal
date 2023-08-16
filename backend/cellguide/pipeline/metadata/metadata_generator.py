from pronto import Ontology

from backend.cellguide.pipeline.metadata.types import CellMetadata
from backend.wmg.data.constants import CL_BASIC_OBO_NAME
from backend.wmg.data.utils import get_pinned_ontology_url


def generate_cellguide_card_metadata(all_cell_type_ids_in_corpus: list[str]) -> dict[str, CellMetadata]:
    ontology = Ontology(get_pinned_ontology_url(CL_BASIC_OBO_NAME))

    cellguide_card_metadata: dict[str, CellMetadata] = {}
    for id in all_cell_type_ids_in_corpus:
        cell = ontology[id]
        if not cell.id.startswith("CL:"):
            continue
        if cell.obsolete:
            continue

        metadata = CellMetadata(
            name=cell.name,
            id=cell.id,
            clDescription=str(cell.definition) if cell.definition is not None else None,
            synonyms=[s.description for s in cell.synonyms],
        )
        cellguide_card_metadata[cell.id] = metadata

    return cellguide_card_metadata
