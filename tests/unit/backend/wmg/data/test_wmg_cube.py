import os.path
import unittest

from backend.wmg.data import tiledb
from backend.wmg.data.schemas.cube_schema import cell_counts_schema
from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME
from backend.wmg.data.wmg_cube import create_cell_count_cube
from tests.unit.backend.wmg.fixtures.test_snapshot import create_temp_wmg_snapshot


@unittest.skip(
    "Requires https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/"
    "chanzuckerberg/single-cell-data-portal/2280"
)
class TestWmgCube(unittest.TestCase):
    @unittest.skip
    def test_create_cube(self):
        # TODO write this test when cube creation code fully implemented
        raise NotImplementedError

    def test__create_cell_counts_cube__sets_schema_correctly(self):
        with create_temp_wmg_snapshot() as snapshot:
            tdb_group = os.path.dirname(snapshot.expression_summary_cube.uri)
            create_cell_count_cube(tdb_group)

        with tiledb.open(f"{tdb_group}/{CELL_COUNTS_CUBE_NAME}") as cell_count_cube:
            self.assertEqual(cell_count_cube.schema, cell_counts_schema)

    @unittest.skip
    def test_create_cube_gracefully_alerts_on_failure(self):
        # TODO write this test when cube creation error handling fully implemented
        raise NotImplementedError
