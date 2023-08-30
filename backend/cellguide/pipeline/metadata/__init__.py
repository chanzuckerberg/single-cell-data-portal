from backend.cellguide.pipeline.constants import CELL_GUIDE_METADATA_FILENAME, CELL_GUIDE_TISSUE_METADATA_FILENAME
from backend.cellguide.pipeline.metadata.metadata_generator import (
    generate_cellguide_card_metadata,
    generate_cellguide_tissue_card_metadata,
)
from backend.cellguide.pipeline.metadata.types import CellMetadata, TissueMetadata
from backend.cellguide.pipeline.ontology_tree import get_ontology_tree_builder
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import output_json
from backend.wmg.api.wmg_api_config import WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.snapshot import load_snapshot


def run(*, output_directory: str):
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    ontology_tree = get_ontology_tree_builder(snapshot=snapshot)
    cell_metadata = get_cell_metadata(ontology_tree=ontology_tree)
    tissue_metadata = get_tissue_metadata(ontology_tree=ontology_tree)

    output_json(cell_metadata, f"{output_directory}/{CELL_GUIDE_METADATA_FILENAME}")
    output_json(tissue_metadata, f"{output_directory}/{CELL_GUIDE_TISSUE_METADATA_FILENAME}")


def get_cell_metadata(*, ontology_tree: OntologyTreeBuilder) -> dict[str, CellMetadata]:
    """
    For all cell type ids in the corpus, this pipeline will generate metadata about each cell, including:
    - name, ex: "native cell"
    - id, ex: "CL:0000003"
    - clDescription, ex: "A cell that is found in a natural setting, which includes multicellular organism cells 'in vivo'
      (i.e. part of an organism), and unicellular organisms 'in environment' (i.e. part of a natural environment)."
    - synonyms, ex: ["cell in vivo"]

    Note that we will be filtering out obsolete cell types and invalid non-CL cell types.
    """
    return generate_cellguide_card_metadata(ontology_tree.all_cell_type_ids_in_corpus)


def get_tissue_metadata(*, ontology_tree: OntologyTreeBuilder) -> dict[str, TissueMetadata]:
    """
    For all tissue ids in the corpus, this pipeline will generate metadata about each tissue, including:
    - name, ex: "lung"
    - id, ex: "UBERON:0002048"
    - uberonDescription, ex: "Respiration organ that develops as an outpocketing of the esophagus."
    - synonyms, ex: ["pulmo"]

    Note that we will be filtering out obsolete tissues.
    """
    return generate_cellguide_tissue_card_metadata(ontology_tree.all_tissue_ids_in_corpus)
