import logging
import warnings

from pronto import Ontology

from backend.cellguide.pipeline.constants import UBERON_BASIC_PERMANENT_URL_PRONTO
from backend.cellguide.pipeline.metadata.types import CellMetadata, TissueMetadata
from backend.wmg.data.constants import CL_BASIC_OBO_NAME
from backend.wmg.data.utils import get_pinned_ontology_url

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
    ontology = Ontology(get_pinned_ontology_url(CL_BASIC_OBO_NAME))

    cellguide_card_metadata: dict[str, CellMetadata] = {}

    obsolete_cell_ids: list[str] = []
    invalid_non_cl_cell_ids: list[str] = []
    cell_ids_with_cl_description = 0
    cell_ids_without_cl_description = 0

    for id in all_cell_type_ids_in_corpus:
        cell = ontology[id]
        if not cell.id.startswith("CL:"):
            invalid_non_cl_cell_ids.append(cell.id)
        elif cell.obsolete:
            obsolete_cell_ids.append(cell.id)
        else:
            description = None
            if cell.definition is not None:
                description = str(cell.definition)
                cell_ids_with_cl_description += 1
            else:
                cell_ids_without_cl_description += 1

            metadata = CellMetadata(
                name=cell.name,
                id=cell.id,
                clDescription=description,
                synonyms=[s.description for s in cell.synonyms],
            )
            cellguide_card_metadata[cell.id] = metadata

    logger.info(f"Filtering out {len(obsolete_cell_ids)} obsolete cell types: {obsolete_cell_ids}")
    logger.info(f"Filtering out {len(invalid_non_cl_cell_ids)} invalid non-CL cell types: {invalid_non_cl_cell_ids}")
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
    with warnings.catch_warnings():
        # loading uberon ontology has some warnings that we don't care about
        warnings.simplefilter("ignore")
        ontology = Ontology(UBERON_BASIC_PERMANENT_URL_PRONTO)

    cellguide_tissue_card_metadata: dict[str, TissueMetadata] = {}

    obsolete_uberon_ids: list[str] = []
    invalid_non_uberon_ids: list[str] = []
    uberon_ids_with_description = 0
    uberon_ids_without_description = 0

    for id in all_tissue_ids_in_corpus:
        tissue = ontology[id]
        if not tissue.id.startswith("UBERON:"):
            invalid_non_uberon_ids.append(tissue.id)
        elif tissue.obsolete:
            obsolete_uberon_ids.append(tissue.id)
        else:
            description = None
            if tissue.definition is not None:
                description = str(tissue.definition)
                uberon_ids_with_description += 1
            else:
                uberon_ids_without_description += 1

            metadata = TissueMetadata(
                name=tissue.name,
                id=tissue.id,
                uberonDescription=description,
                synonyms=[s.description for s in tissue.synonyms],
            )
            cellguide_tissue_card_metadata[tissue.id] = metadata

    logger.info(f"Filtering out {len(obsolete_uberon_ids)} obsolete tissues: {obsolete_uberon_ids}")
    logger.info(f"Filtering out {len(invalid_non_uberon_ids)} invalid non-UBERON tissues: {invalid_non_uberon_ids}")
    logger.info(
        f"Found {uberon_ids_with_description} tissues with UBERON descriptions and {uberon_ids_without_description} without UBERON descriptions"
    )

    return cellguide_tissue_card_metadata
