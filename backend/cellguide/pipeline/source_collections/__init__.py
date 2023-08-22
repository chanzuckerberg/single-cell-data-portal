from backend.cellguide.pipeline.constants import SOURCE_COLLECTIONS_FOLDERNAME
from backend.cellguide.pipeline.ontology_tree import get_ontology_tree_builder
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.source_collections.source_collections_generator import generate_source_collections_data
from backend.cellguide.pipeline.source_collections.types import SourceCollectionsData
from backend.cellguide.pipeline.utils import output_json_per_key
from backend.wmg.api.wmg_api_config import WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.snapshot import load_snapshot


def run(*, output_directory: str) -> dict[str, SourceCollectionsData]:
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    ontology_tree = get_ontology_tree_builder(snapshot=snapshot)
    data = get_source_collections_data(ontology_tree=ontology_tree)
    output_json_per_key(data, f"{output_directory}/{SOURCE_COLLECTIONS_FOLDERNAME}")


def get_source_collections_data(*, ontology_tree: OntologyTreeBuilder):
    # this one-line wrapper is used to match the pattern used in the other sub-pipelines
    # for convenient imports
    return generate_source_collections_data(ontology_tree.all_cell_type_ids_in_corpus)
