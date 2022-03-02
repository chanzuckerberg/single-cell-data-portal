import json
import unittest
from unittest.mock import patch

from backend.corpora.api_server.app import app
from backend.wmg.api import v1
from backend.wmg.data.schema import cube_non_indexed_dims
from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.wmg.fixtures.cube import create_temp_cube, all_ones_attr_values


class WmgApiV1Tests(unittest.TestCase):
    """
    Tests WMG API endpoints. Tests the flask app only, and not other stack dependencies, such as S3. Builds and uses a
    temporary WMG cube on local filesystem to avoid dependency on localstack S3.
    """

    def setUp(self):
        super().setUp()
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            self.app = app.test_client(use_cookies=False)

    def test__primary_filter_dimensions__returns_200(self):
        response = self.app.get("/wmg/v1/primary_filter_dimensions")

        self.assertEqual(200, response.status_code)

    def test__primary_filter_dimensions__returns_valid_response_body(self):
        response = self.app.get("/wmg/v1/primary_filter_dimensions")

        expected = dict(
            snapshot_id=v1.DUMMY_SNAPSHOT_UUID,
            organism_terms=[dict(oid1="olbl1"), dict(oid2="olbl2")],
            tissue_terms=[dict(ttid1="ttlbl1"), dict(ttid2="ttlbl2")],
        )

        self.assertEqual(expected, json.loads(response.data))

    @patch("backend.wmg.api.v1.find_cube_latest_snapshot")
    def test__query__minimal_valid_request_returns_200_and_empty_expr_summary(self, find_cube_latest_snapshot):
        with create_temp_cube() as cube:
            # setup up API endpoints to use a cube containing all stat values of 1, for a deterministic expected query
            # response
            find_cube_latest_snapshot.return_value = cube

            request = dict(
                filter=dict(
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                ),
                response_option="include_filter_dims_include_dataset_links",
            )

            response = self.app.post("/wmg/v1/query", json=request)

            self.assertEqual(200, response.status_code)

            expected_response = {
                "snapshot_id": v1.DUMMY_SNAPSHOT_UUID,
                "expression_summary": {},
                "term_id_labels": {"cell_types": [], "genes": []},
            }

            self.assertEqual(expected_response, json.loads(response.data))

    def test__query__empty_request_returns_400(self):
        response = self.app.post("/wmg/v1/query", json={})

        self.assertEqual(400, response.status_code)
        self.assertEqual("'filter' is a required property", json.loads(response.data)["detail"])

    def test__query__missing_organism_request_returns_400(self):
        request = dict(
            filter=dict(tissue_ontology_term_ids=["tissue_ontology_term_id_0"]),
            response_option="include_filter_dims_include_dataset_links",
        )

        response = self.app.post("/wmg/v1/query", json=request)

        self.assertEqual(400, response.status_code)

    def test__query__missing_tissue_request_returns_400(self):
        request = dict(
            filter=dict(organism_ontology_term_id="organism_ontology_term_id_0"),
            response_option="include_filter_dims_include_dataset_links",
        )

        response = self.app.post("/wmg/v1/query", json=request)

        self.assertEqual(400, response.status_code)

    @patch("backend.wmg.api.v1.find_cube_latest_snapshot")
    def test__query__valid_request_returns_valid_response_body(self, find_cube_latest_snapshot):
        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as all_ones_cube:
            # setup up API endpoints to use a cube containing all stat values of 1, for a deterministic expected query
            # response
            find_cube_latest_snapshot.return_value = all_ones_cube

            request = dict(
                filter=dict(
                    gene_ontology_term_ids=[
                        "gene_ontology_term_id_0",
                        "gene_ontology_term_id_2",
                    ],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_1", "tissue_ontology_term_id_2"],
                ),
                response_option="include_filter_dims_include_dataset_links",
            )

            response = self.app.post("/wmg/v1/query", json=request)

            # sanity check the expected values of the `n` stat for each data viz point; if this fails, the cube test
            # fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture has
            # changed
            expected_cell_count_per_cell_type = dim_size ** len(
                set(cube_non_indexed_dims).difference({"cell_type_ontology_term_id"})
            )
            assert expected_cell_count_per_cell_type == 729

            expected = {
                "snapshot_id": v1.DUMMY_SNAPSHOT_UUID,
                "expression_summary": {
                    "gene_ontology_term_id_0": {
                        "tissue_ontology_term_id_1": [
                            {"id": "cell_type_ontology_term_id_0", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                            {"id": "cell_type_ontology_term_id_1", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                            {"id": "cell_type_ontology_term_id_2", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                        ],
                        "tissue_ontology_term_id_2": [
                            {"id": "cell_type_ontology_term_id_0", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                            {"id": "cell_type_ontology_term_id_1", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                            {"id": "cell_type_ontology_term_id_2", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                        ],
                    },
                    "gene_ontology_term_id_2": {
                        "tissue_ontology_term_id_1": [
                            {"id": "cell_type_ontology_term_id_0", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                            {"id": "cell_type_ontology_term_id_1", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                            {"id": "cell_type_ontology_term_id_2", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                        ],
                        "tissue_ontology_term_id_2": [
                            {"id": "cell_type_ontology_term_id_0", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                            {"id": "cell_type_ontology_term_id_1", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                            {"id": "cell_type_ontology_term_id_2", "n": 729, "me": 1.0, "pc": 0.0, "tpc": 0.0},
                        ],
                    },
                },
                "term_id_labels": {
                    "cell_types": [
                        {"cell_type_ontology_term_id_0": "cell_type_ontology_term_id_0_label"},
                        {"cell_type_ontology_term_id_1": "cell_type_ontology_term_id_1_label"},
                        {"cell_type_ontology_term_id_2": "cell_type_ontology_term_id_2_label"},
                    ],
                    "genes": [
                        {"gene_ontology_term_id_0": "gene_ontology_term_id_0_label"},
                        {"gene_ontology_term_id_2": "gene_ontology_term_id_2_label"},
                    ],
                },
            }
            self.assertEqual(expected, json.loads(response.data))
