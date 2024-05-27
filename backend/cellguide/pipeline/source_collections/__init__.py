from backend.cellguide.common.constants import SOURCE_COLLECTIONS_FOLDERNAME
from backend.cellguide.pipeline.source_collections.source_collections_generator import generate_source_collections_data
from backend.cellguide.pipeline.source_collections.types import SourceCollectionsData
from backend.cellguide.pipeline.utils import output_json_per_key
from backend.common.census_cube.utils import get_all_cell_type_ids_in_corpus


def run(*, output_directory: str) -> dict[str, SourceCollectionsData]:
    data = get_source_collections_data()
    output_json_per_key(data, f"{output_directory}/{SOURCE_COLLECTIONS_FOLDERNAME}")


def get_source_collections_data():
    # this one-line wrapper is used to match the pattern used in the other sub-pipelines
    # for convenient imports
    all_cell_type_ids_in_corpus = get_all_cell_type_ids_in_corpus()
    return generate_source_collections_data(all_cell_type_ids_in_corpus)
