import logging

from backend.cellguide.pipeline.metadata.types import CellMetadata, TissueMetadata
from backend.common.census_cube.utils import ontology_parser

logger = logging.getLogger(__name__)


def generate_cellguide_card_metadata(all_cell_type_ids_in_corpus: list[str]) -> dict[str, CellMetadata]:
    """
    For all cell type ids in the corpus, this pipeline will generate metadata about each cell, including:
    - name, ex: "native cell"
    - id, ex: "CL:0000003"
    - clDescription, ex: "A cell that is found in a natural setting, which includes multicellular organism cells 'in vivo'
      (i.e. part of an organism), and unicellular organisms 'in environment' (i.e. part of a natural environment)."
    - synonyms, ex: ["cell in vivo"]

    Note that we will be filtering out obsolete cell types and invalid non-CL cell types.
    """
    logger.info(f"Generating cellguide card metadata for {len(all_cell_type_ids_in_corpus)} cell types...")

    cellguide_card_metadata: dict[str, CellMetadata] = {}

    obsolete_cell_ids: list[str] = []
    cell_ids_with_cl_description = 0
    cell_ids_without_cl_description = 0

    for id in all_cell_type_ids_in_corpus:

        if ontology_parser.is_term_deprecated(id):
            obsolete_cell_ids.append(id)
        else:
            description = ontology_parser.get_term_description(id)
            if description is not None:
                cell_ids_with_cl_description += 1
            else:
                cell_ids_without_cl_description += 1

            metadata = CellMetadata(
                name=ontology_parser.get_term_label(id),
                id=id,
                clDescription=description,
                synonyms=ontology_parser.get_term_synonyms(id),
            )
            cellguide_card_metadata[id] = metadata

    logger.info(f"Filtering out {len(obsolete_cell_ids)} obsolete cell types: {obsolete_cell_ids}")
    logger.info(
        f"Found {cell_ids_with_cl_description} cell types with CL descriptions and {cell_ids_without_cl_description} without CL descriptions"
    )

    return cellguide_card_metadata


def generate_cellguide_tissue_card_metadata(all_tissue_ids_in_corpus: list[str]) -> dict[str, CellMetadata]:
    """
    For all tissue ids in the corpus, this pipeline will generate metadata about each tissue, including:
    - name, ex: "lung"
    - id, ex: "UBERON:0002048"
    - uberonDescription, ex: "Respiration organ that develops as an outpocketing of the esophagus."
    - synonyms, ex: ["pulmo"]

    Note that we will be filtering out obsolete tissues.
    """
    logger.info(f"Generating cellguide tissue card metadata for {len(all_tissue_ids_in_corpus)} tissues...")

    cellguide_tissue_card_metadata: dict[str, TissueMetadata] = {}

    obsolete_uberon_ids: list[str] = []
    uberon_ids_with_description = 0
    uberon_ids_without_description = 0

    for id in all_tissue_ids_in_corpus:
        if ontology_parser.is_term_deprecated(id):
            obsolete_uberon_ids.append(id)
        else:
            description = ontology_parser.get_term_description(id)
            if description is not None:
                uberon_ids_with_description += 1
            else:
                uberon_ids_without_description += 1

            metadata = TissueMetadata(
                name=ontology_parser.get_term_label(id),
                id=id,
                uberonDescription=description,
                synonyms=ontology_parser.get_term_synonyms(id),
            )
            cellguide_tissue_card_metadata[id] = metadata

    logger.info(f"Filtering out {len(obsolete_uberon_ids)} obsolete tissues: {obsolete_uberon_ids}")
    logger.info(
        f"Found {uberon_ids_with_description} tissues with UBERON descriptions and {uberon_ids_without_description} without UBERON descriptions"
    )

    return cellguide_tissue_card_metadata
