from backend.cellguide.common.constants import SOURCE_COLLECTIONS_FOLDERNAME
from backend.cellguide.pipeline.constants import CELLGUIDE_CENSUS_CUBE_DATA_SCHEMA_VERSION
from backend.cellguide.pipeline.source_collections.source_collections_generator import generate_source_collections_data
from backend.cellguide.pipeline.source_collections.types import SourceCollectionsData
from backend.cellguide.pipeline.utils import output_json_per_key
from backend.common.census_cube.data import snapshot as sn
from backend.common.census_cube.utils import get_all_cell_type_ids_in_corpus


def run(*, output_directory: str) -> dict[str, SourceCollectionsData]:
    data = get_source_collections_data()
    output_json_per_key(data, f"{output_directory}/{SOURCE_COLLECTIONS_FOLDERNAME}")


def get_source_collections_data():
    # this one-line wrapper is used to match the pattern used in the other sub-pipelines
    # for convenient imports
    snapshot = sn.load_snapshot(snapshot_schema_version=CELLGUIDE_CENSUS_CUBE_DATA_SCHEMA_VERSION)
    all_cell_type_ids_in_corpus = get_all_cell_type_ids_in_corpus(snapshot)
    return generate_source_collections_data(all_cell_type_ids_in_corpus)
