import os
import sys
from tempfile import TemporaryDirectory

from backend.cellguide.pipeline.canonical_marker_genes import get_canonical_marker_genes
from backend.cellguide.pipeline.computational_marker_genes import get_computational_marker_genes
from backend.cellguide.pipeline.metadata import get_cell_metadata, get_tissue_metadata
from backend.cellguide.pipeline.ontology_tree import get_ontology_tree_data
from backend.cellguide.pipeline.source_collections import get_source_collections_data
from backend.cellguide.pipeline.utils import output_json
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot
from tests.unit.cellguide_pipeline.constants import (
    CANONICAL_MARKER_GENES_FIXTURE_FILENAME,
    CELLGUIDE_PIPELINE_FIXTURES_BASEPATH,
    CELLTYPE_METADATA_FIXTURE_FILENAME,
    CELLTYPE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME,
    COMPUTATIONAL_MARKER_GENES_FIXTURE_FILENAME,
    ONTOLOGY_GRAPH_FIXTURE_FILENAME,
    SOURCE_COLLECTIONS_FIXTURE_FILENAME,
    TISSUE_METADATA_FIXTURE_FILENAME,
    TISSUE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME,
)

root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)

TEST_SNAPSHOT = "realistic-test-snapshot"


""" ################################# DANGER #######################################

Run this file if and only if you are confident there are no bugs in the CellGuide
pipeline.

Any and all unit test assertion errors prior to running this script must be expected
due to intended changes in the pipeline.

------------------------------------------------------------------------------------

This module generates the CellGuide data using the test snapshot stored in {TEST_SNAPSHOT}.
Requires an internet connection.

Run this script with the below command from the root directory of this repo:
```
python -m scripts.generate_cellguide_pipeline_test_fixtures
```
"""


def run_cellguide_pipeline():
    with TemporaryDirectory(), TemporaryDirectory(), load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
        # Get ontology tree data
        ontology_tree_data = get_ontology_tree_data(snapshot=snapshot)
        ontology_tree = ontology_tree_data.tree_builder

        # Get cell metadata
        cell_metadata = get_cell_metadata(ontology_tree=ontology_tree)

        # Get tissue metadata
        tissue_metadata = get_tissue_metadata(ontology_tree=ontology_tree)

        # Get canonical marker genes
        canonical_marker_genes = get_canonical_marker_genes(snapshot=snapshot, ontology_tree=ontology_tree)

        # Get source data
        source_collections = get_source_collections_data(ontology_tree=ontology_tree)

        # Get computatoinal marker genes
        computational_marker_genes = get_computational_marker_genes(snapshot=snapshot, ontology_tree=ontology_tree)

        output_json(
            ontology_tree_data.ontology_graph,
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ONTOLOGY_GRAPH_FIXTURE_FILENAME}",
        )
        output_json(
            ontology_tree_data.all_states_per_cell_type,
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{CELLTYPE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME}",
        )
        output_json(
            ontology_tree_data.all_states_per_tissue,
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{TISSUE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME}",
        )
        output_json(
            cell_metadata,
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{CELLTYPE_METADATA_FIXTURE_FILENAME}",
        )
        output_json(
            tissue_metadata,
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{TISSUE_METADATA_FIXTURE_FILENAME}",
        )
        output_json(
            canonical_marker_genes,
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{CANONICAL_MARKER_GENES_FIXTURE_FILENAME}",
        )
        output_json(
            computational_marker_genes,
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{COMPUTATIONAL_MARKER_GENES_FIXTURE_FILENAME}",
        )

        output_json(
            source_collections,
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{SOURCE_COLLECTIONS_FIXTURE_FILENAME}",
        )


if __name__ == "__main__":
    run_cellguide_pipeline()
