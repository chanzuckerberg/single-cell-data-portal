import logging

from pronto import Ontology

from backend.geneguide.pipeline.constants import GO_URL
from backend.geneguide.pipeline.metadata.types import GoTermMetadata

logger = logging.getLogger(__name__)


def generate_card_metadata(all_go_term_ids: list[str]) -> dict[str, GoTermMetadata]:
    """
    For all GO terms, this pipeline will generate metadata about each, including:
    - name: str,
    - id: str,
    - goDescription: Optional[str],
    - synonyms: list[str],

    Note that we will be filtering out obsolete GO terms.
    """
    logger.info(f"Generating card metadata for {len(all_go_term_ids)} GO terms...")
    ontology = Ontology(GO_URL)

    card_metadata: dict[str, GoTermMetadata] = {}

    obsolete_go_terms: list[str] = []
    go_terms_with_description = 0
    go_terms_without_description = 0

    for id in all_go_term_ids:
        goterm = ontology[id]
        if not goterm.id.startswith("GO:"):
            continue

        if goterm.obsolete:
            obsolete_go_terms.append(goterm.id)
        else:
            description = None
            if goterm.definition is not None:
                description = str(goterm.definition)
                go_terms_with_description += 1
            else:
                go_terms_without_description += 1

            metadata = GoTermMetadata(
                name=goterm.name,
                id=goterm.id,
                goDescription=description,
                synonyms=[s.description for s in goterm.synonyms],
            )
            card_metadata[goterm.id] = metadata

    logger.info(f"Filtering out {len(obsolete_go_terms)} obsolete cell types: {obsolete_go_terms}")
    logger.info(
        f"Found {go_terms_with_description} cell types with GO descriptions and {go_terms_without_description} without GO descriptions"
    )

    return card_metadata
