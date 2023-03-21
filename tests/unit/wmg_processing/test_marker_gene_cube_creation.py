import os
import shutil
import unittest

import tiledb
from pandas.testing import assert_series_equal

from backend.wmg.pipeline.summary_cubes.marker_genes import create_marker_genes_cube
from tests.unit.backend.wmg.fixtures import FIXTURES_ROOT

TEST_SNAPSHOT = "realistic-test-snapshot"

# creates the marker gene cube de novo and compares to the existing fixture
# note, this test only implicitly tests the marker gene computation logic
# by ensuring that the generated cube matches that provided in the snapshot fixture
# the calculation logic is being tested directly in test_marker_gene_calculation.py

# to regenerate the test FMG snapshot fixture from an existing (full) snapshot, please run:
# python -m scripts.generate_fmg_test_snapshot "path/to/full/snapshot" "path/to/new/fixture"


class MarkerGeneCubeCreationTest(unittest.TestCase):
    def setUp(self):
        super().setUp()
        os.rename(f"{FIXTURES_ROOT}/{TEST_SNAPSHOT}/marker_genes", f"{FIXTURES_ROOT}/{TEST_SNAPSHOT}/marker_genes_old")

    def test__marker_gene_cube_creation(self):
        create_marker_genes_cube(f"{FIXTURES_ROOT}/{TEST_SNAPSHOT}")
        with tiledb.open(f"{FIXTURES_ROOT}/{TEST_SNAPSHOT}/marker_genes") as new_cube, tiledb.open(
            f"{FIXTURES_ROOT}/{TEST_SNAPSHOT}/marker_genes_old"
        ) as ref_cube:
            self.assertEqual(str(new_cube.schema), str(ref_cube.schema))
            df1 = (
                ref_cube.df[:]
                .sort_values(["cell_type_ontology_term_id", "gene_ontology_term_id"])
                .reset_index(drop=True)
            )
            df2 = (
                new_cube.df[:]
                .sort_values(["cell_type_ontology_term_id", "gene_ontology_term_id"])
                .reset_index(drop=True)
            )
            [assert_series_equal(df1[col], df2[col], rtol=1e-3) for col in df1]

    def tearDown(self):
        if os.path.exists(f"{FIXTURES_ROOT}/{TEST_SNAPSHOT}/marker_genes"):
            shutil.rmtree(f"{FIXTURES_ROOT}/{TEST_SNAPSHOT}/marker_genes")
        os.rename(f"{FIXTURES_ROOT}/{TEST_SNAPSHOT}/marker_genes_old", f"{FIXTURES_ROOT}/{TEST_SNAPSHOT}/marker_genes")
        super().tearDown()
