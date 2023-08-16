from backend.cellguide.pipeline.constants import CELL_GUIDE_METADATA_FILENAME
from backend.cellguide.pipeline.metadata.metadata_generator import generate_cellguide_card_metadata
from backend.cellguide.pipeline.metadata.types import CellMetadata
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import output_json


def run(*, output_directory: str, ontology_tree: OntologyTreeBuilder) -> dict[str, CellMetadata]:
    """
    For all cell type ids in the corpus, this pipeline will generate metadata about each cell, including:
    - name, ex: "native cell"
    - id, ex: "CL:0000003"
    - clDescription, ex: "A cell that is found in a natural setting, which includes multicellular organism cells 'in vivo'
      (i.e. part of an organism), and unicellular organisms 'in environment' (i.e. part of a natural environment)."
    - synonyms, ex: ["cell in vivo"]

    Note that we will be filtering out obsolete cell types and invalid non-CL cell types.
    """
    cell_metadata = generate_cellguide_card_metadata(ontology_tree.all_cell_type_ids_in_corpus)
    output_json(cell_metadata, f"{output_directory}/{CELL_GUIDE_METADATA_FILENAME}")
    return cell_metadata
