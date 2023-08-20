import os
import shutil
import sys
import tempfile
import unittest.mock
from functools import partial

from backend.cellguide.pipeline.canonical_marker_genes import run as run_canonical_marker_gene_pipeline
from backend.cellguide.pipeline.computational_marker_genes import (
    get_marker_genes_per_and_across_tissues as get_computational_marker_genes,
)
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

root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)

""" ################################# DANGER #################################### 

Run this file if and only if you are confident there are no bugs in the CellGuide
pipeline.

Any and all unit test assertion errors prior to running this script must be expected
due to intended changes in the pipeline.
"""

TEST_SNAPSHOT = "realistic-test-snapshot"

PATH_TO_FIXTURES = "tests/unit/cellguide_pipeline/fixtures/"


def custom_load_snapshot_into_cube_dir(cube_dir: str, **kwargs):
    return load_realistic_test_snapshot_obj(TEST_SNAPSHOT, cube_dir)


def run_cellguide_pipeline():
    output_directory = sys.argv[1]
    with tempfile.TemporaryDirectory() as cube_dir:
        # Patch load_snapshot with custom_load_snapshot
        custom_load_snapshot = partial(custom_load_snapshot_into_cube_dir, cube_dir=cube_dir)
        with unittest.mock.patch("backend.cellguide.pipeline.ontology_tree.load_snapshot", new=custom_load_snapshot):
            # Run ontology tree pipeline
            ontology_tree = run_ontology_tree_pipeline(output_directory)

        with unittest.mock.patch("backend.cellguide.pipeline.metadata.load_snapshot", new=custom_load_snapshot):
            # Generate cell guide cards, synonyms, and descriptions
            run_metadata_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

        with unittest.mock.patch(
            "backend.cellguide.pipeline.canonical_marker_genes.load_snapshot", new=custom_load_snapshot
        ):
            # Generate canonical marker genes from ASCT-B (HUBMAP)
            run_canonical_marker_gene_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

        with unittest.mock.patch(
            "backend.cellguide.pipeline.source_collections.load_snapshot", new=custom_load_snapshot
        ):
            # Generate source data for each cell type
            run_source_collections_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)

        with unittest.mock.patch(
            "backend.cellguide.pipeline.computational_marker_genes.load_snapshot", new=custom_load_snapshot
        ):
            # Generate computational marker genes from the CZI corpus
            marker_genes = get_computational_marker_genes(
                output_directory=output_directory, ontology_tree=ontology_tree
            )
            output_json(marker_genes, f"{output_directory}/computational_marker_genes.json")

    shutil.move(f"{output_directory}/{ONTOLOGY_TREE_FILENAME}", f"{output_directory}/ontology_graph.json")
    shutil.move(
        f"{output_directory}/{ONTOLOGY_TREE_STATE_PER_CELLTYPE_FILENAME}",
        f"{output_directory}/all_states_per_cell_type.json",
    )
    shutil.move(
        f"{output_directory}/{ONTOLOGY_TREE_STATE_PER_TISSUE_FILENAME}",
        f"{output_directory}/all_states_per_tissue.json",
    )
    shutil.move(f"{output_directory}/{SOURCE_COLLECTIONS_FILENAME}", f"{output_directory}/source_collections.json")
    shutil.move(f"{output_directory}/{CELL_GUIDE_METADATA_FILENAME}", f"{output_directory}/cell_metadata.json")
    shutil.move(f"{output_directory}/{CELL_GUIDE_TISSUE_METADATA_FILENAME}", f"{output_directory}/tissue_metadata.json")
    shutil.move(
        f"{output_directory}/{CANONICAL_MARKER_GENES_FILENAME}", f"{output_directory}/canonical_marker_genes.json"
    )
    shutil.move(output_directory, PATH_TO_FIXTURES)


if __name__ == "__main__":
    run_cellguide_pipeline()
