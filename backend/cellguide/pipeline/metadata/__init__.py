from backend.cellguide.common.constants import CELL_GUIDE_METADATA_FILENAME, CELL_GUIDE_TISSUE_METADATA_FILENAME
from backend.cellguide.pipeline.constants import CELLGUIDE_CENSUS_CUBE_DATA_SCHEMA_VERSION
from backend.cellguide.pipeline.metadata.metadata_generator import (
    generate_cellguide_card_metadata,
    generate_cellguide_tissue_card_metadata,
)
from backend.cellguide.pipeline.metadata.types import CellMetadata, TissueMetadata
from backend.cellguide.pipeline.utils import output_json
from backend.common.census_cube.data import snapshot as sn
from backend.common.census_cube.utils import get_all_cell_type_ids_in_corpus, get_all_tissue_ids_in_corpus


def run(*, output_directory: str):
    cell_metadata = get_cell_metadata()
    tissue_metadata = get_tissue_metadata()

    output_json(cell_metadata, f"{output_directory}/{CELL_GUIDE_METADATA_FILENAME}")
    output_json(tissue_metadata, f"{output_directory}/{CELL_GUIDE_TISSUE_METADATA_FILENAME}")


def get_cell_metadata() -> dict[str, CellMetadata]:
    """
    For all cell type ids in the corpus, this pipeline will generate metadata about each cell, including:
    - name, ex: "native cell"
    - id, ex: "CL:0000003"
    - clDescription, ex: "A cell that is found in a natural setting, which includes multicellular organism cells 'in vivo'
      (i.e. part of an organism), and unicellular organisms 'in environment' (i.e. part of a natural environment)."
    - synonyms, ex: ["cell in vivo"]

    Note that we will be filtering out obsolete cell types and invalid non-CL cell types.
    """
    snapshot = sn.load_snapshot(snapshot_schema_version=CELLGUIDE_CENSUS_CUBE_DATA_SCHEMA_VERSION)
    all_cell_type_ids_in_corpus = get_all_cell_type_ids_in_corpus(snapshot)
    return generate_cellguide_card_metadata(all_cell_type_ids_in_corpus)


def get_tissue_metadata() -> dict[str, TissueMetadata]:
    """
    For all tissue ids in the corpus, this pipeline will generate metadata about each tissue, including:
    - name, ex: "lung"
    - id, ex: "UBERON:0002048"
    - uberonDescription, ex: "Respiration organ that develops as an outpocketing of the esophagus."
    - synonyms, ex: ["pulmo"]

    Note that we will be filtering out obsolete tissues.
    """
    snapshot = sn.load_snapshot(snapshot_schema_version=CELLGUIDE_CENSUS_CUBE_DATA_SCHEMA_VERSION)
    all_tissue_ids_in_corpus = get_all_tissue_ids_in_corpus(snapshot)
    return generate_cellguide_tissue_card_metadata(all_tissue_ids_in_corpus)
