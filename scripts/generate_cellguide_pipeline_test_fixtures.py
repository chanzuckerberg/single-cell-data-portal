import os
import shutil
import sys
import unittest.mock
from functools import partial
from tempfile import TemporaryDirectory

from backend.cellguide.pipeline.canonical_marker_genes import run as run_canonical_marker_gene_pipeline
from backend.cellguide.pipeline.computational_marker_genes import run as run_computational_marker_genes_pipeline
from backend.cellguide.pipeline.computational_marker_genes.types import ComputationalMarkerGenes
from backend.cellguide.pipeline.constants import (
    CANONICAL_MARKER_GENES_FILENAME,
    CELL_GUIDE_METADATA_FILENAME,
    CELL_GUIDE_TISSUE_METADATA_FILENAME,
    ONTOLOGY_TREE_FILENAME,
    ONTOLOGY_TREE_STATE_PER_CELLTYPE_FILENAME,
    ONTOLOGY_TREE_STATE_PER_TISSUE_FILENAME,
    SOURCE_COLLECTIONS_FILENAME,
)
from backend.cellguide.pipeline.metadata import run as run_metadata_pipeline
from backend.cellguide.pipeline.ontology_tree import run as run_ontology_tree_pipeline
from backend.cellguide.pipeline.source_collections import run as run_source_collections_pipeline
from backend.cellguide.pipeline.utils import output_json
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot_obj
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


def custom_load_snapshot_into_cube_dir(cube_dir: str, **kwargs):
    return load_realistic_test_snapshot_obj(TEST_SNAPSHOT, cube_dir)


def custom_output_computational_marker_genes(
    marker_genes: dict[str, list[ComputationalMarkerGenes]], output_directory: str
):
    output_json(marker_genes, f"{output_directory}/{COMPUTATIONAL_MARKER_GENES_FIXTURE_FILENAME}")


def run_cellguide_pipeline():
    with TemporaryDirectory() as cube_dir, TemporaryDirectory() as output_directory:
        # Patch load_snapshot with custom_load_snapshot
        custom_load_snapshot = partial(custom_load_snapshot_into_cube_dir, cube_dir=cube_dir)
        with unittest.mock.patch("backend.cellguide.pipeline.ontology_tree.load_snapshot", new=custom_load_snapshot):
            # Run ontology tree pipeline
            ontology_tree = run_ontology_tree_pipeline(output_directory)

        run_metadata_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

        with unittest.mock.patch(
            "backend.cellguide.pipeline.canonical_marker_genes.load_snapshot", new=custom_load_snapshot
        ):
            # Generate canonical marker genes from ASCT-B (HUBMAP)
            run_canonical_marker_gene_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

        # Generate source data for each cell type
        run_source_collections_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

        with unittest.mock.patch(
            "backend.cellguide.pipeline.computational_marker_genes.load_snapshot", new=custom_load_snapshot
        ), unittest.mock.patch(
            "backend.cellguide.pipeline.computational_marker_genes.output_marker_genes",
            new=custom_output_computational_marker_genes,
        ):
            # Generate computational marker genes from the CZI corpus
            run_computational_marker_genes_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

        shutil.move(
            f"{output_directory}/{ONTOLOGY_TREE_FILENAME}",
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{ONTOLOGY_GRAPH_FIXTURE_FILENAME}",
        )
        shutil.move(
            f"{output_directory}/{ONTOLOGY_TREE_STATE_PER_CELLTYPE_FILENAME}",
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{CELLTYPE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME}",
        )
        shutil.move(
            f"{output_directory}/{ONTOLOGY_TREE_STATE_PER_TISSUE_FILENAME}",
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{TISSUE_ONTOLOGY_TREE_STATE_FIXTURE_FILENAME}",
        )
        shutil.move(
            f"{output_directory}/{SOURCE_COLLECTIONS_FILENAME}",
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{SOURCE_COLLECTIONS_FIXTURE_FILENAME}",
        )
        shutil.move(
            f"{output_directory}/{CELL_GUIDE_METADATA_FILENAME}",
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{CELLTYPE_METADATA_FIXTURE_FILENAME}",
        )
        shutil.move(
            f"{output_directory}/{CELL_GUIDE_TISSUE_METADATA_FILENAME}",
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{TISSUE_METADATA_FIXTURE_FILENAME}",
        )
        shutil.move(
            f"{output_directory}/{CANONICAL_MARKER_GENES_FILENAME}",
            f"{CELLGUIDE_PIPELINE_FIXTURES_BASEPATH}/{CANONICAL_MARKER_GENES_FIXTURE_FILENAME}",
        )


if __name__ == "__main__":
    run_cellguide_pipeline()
