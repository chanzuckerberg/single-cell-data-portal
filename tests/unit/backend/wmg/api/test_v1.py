import json
import unittest
from unittest.mock import patch

from backend.api_server.app import app
from backend.wmg.data.schemas.cube_schema import expression_summary_non_indexed_dims
from tests.unit.backend.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.wmg.fixtures.test_primary_filters import (
    test_gene_terms,
    test_organism_terms,
    test_snapshot_id,
    test_tissue_terms,
)
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    all_ones_expression_summary_values,
    all_tens_cell_counts_values,
    all_X_cell_counts_values,
    create_temp_wmg_snapshot,
    exclude_all_but_one_gene_per_organism,
    load_realistic_test_snapshot,
    reverse_cell_type_ordering,
)

TEST_SNAPSHOT = "realistic-test-snapshot"


class WmgApiV1Tests(unittest.TestCase):
    """
    Tests WMG API endpoints. Tests the flask app only, and not other stack dependencies, such as S3. Builds and uses a
    temporary WMG cube on local filesystem to avoid dependency on localstack S3.
    """

    def setUp(self):
        super().setUp()
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            self.app = app.test_client(use_cookies=False)

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.maxDiff = None

    @patch("backend.wmg.api.v1.load_snapshot")
    def test__primary_filter_dimensions__returns_200(self, load_snapshot):
        # This test appears to be hitting a TileDB (<=0.13.1) bug and fails (intermittently) if dim_size=3
        with create_temp_wmg_snapshot(dim_size=1) as snapshot:
            # setup up API endpoints to use a mocked cube
            load_snapshot.return_value = snapshot
            response = self.app.get("/wmg/v1/primary_filter_dimensions")

        self.assertEqual(200, response.status_code)

    @patch("backend.wmg.api.v1.gene_term_label")
    @patch("backend.wmg.api.v1.ontology_term_label")
    @patch("backend.wmg.api.v1.load_snapshot")
    def test__primary_filter_dimensions__returns_valid_response_body(
        self, load_snapshot, ontology_term_label, gene_term_label
    ):
        with create_temp_wmg_snapshot(
            dim_size=3, exclude_logical_coord_fn=exclude_all_but_one_gene_per_organism
        ) as snapshot:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_snapshot.return_value = snapshot

            # mock the functions in the ontology_labels module, so we can assert deterministic values in the
            # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
            # module is separately unit tested, and here we just want to verify the response building logic is correct.
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"

            response = self.app.get("/wmg/v1/primary_filter_dimensions")

        expected = dict(
            snapshot_id=test_snapshot_id,
            organism_terms=test_organism_terms,
            tissue_terms=test_tissue_terms,
            gene_terms=test_gene_terms,
        )

        self.assertEqual(expected, json.loads(response.data))

    @patch("backend.wmg.api.v1.gene_term_label")
    @patch("backend.wmg.api.v1.ontology_term_label")
    @patch("backend.wmg.api.v1.load_snapshot")
    def test__query_single_primary_dims__returns_200_and_correct_response(
        self, load_snapshot, ontology_term_label, gene_term_label
    ):
        dim_size = 1
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_snapshot.return_value = snapshot

            # mock the functions in the ontology_labels module, so we can assert deterministic values in the
            # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
            # module is separately unit tested, and here we just want to verify the response building logic is correct.
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"

            request = dict(
                filter=dict(
                    gene_ontology_term_ids=[
                        "gene_ontology_term_id_0",
                    ],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                ),
            )

            response = self.app.post("/wmg/v1/query", json=request)

            self.assertEqual(200, response.status_code)

            expected_response = {
                "snapshot_id": "dummy-snapshot",
                "expression_summary": {
                    "gene_ontology_term_id_0": {
                        "tissue_ontology_term_id_0": [
                            {"id": "cell_type_ontology_term_id_0", "me": 1.0, "n": 1, "pc": 0.1, "tpc": 0.1}
                        ]
                    },
                },
                "term_id_labels": {
                    "cell_types": {
                        "tissue_ontology_term_id_0": [
                            {
                                "cell_type": "cell_type_ontology_term_id_0_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                                "total_count": 10,
                                "depth": 0,
                            }
                        ]
                    },
                    "genes": [{"gene_ontology_term_id_0": "gene_ontology_term_id_0_label"}],
                },
                "filter_dims": {},
            }
            self.assertEqual(expected_response, json.loads(response.data))

    @patch("backend.wmg.api.v1.gene_term_label")
    @patch("backend.wmg.api.v1.ontology_term_label")
    @patch("backend.wmg.api.v1.load_snapshot")
    def test__query_request_multi_primary_dims_only__returns_200_and_correct_response(
        self, load_snapshot, ontology_term_label, gene_term_label
    ):
        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_snapshot.return_value = snapshot

            # mock the functions in the ontology_labels module, so we can assert deterministic values in the
            # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
            # module is separately unit tested, and here we just want to verify the response building logic is correct.
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"

            request = dict(
                filter=dict(
                    gene_ontology_term_ids=[
                        "gene_ontology_term_id_0",
                        "gene_ontology_term_id_2",
                    ],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_1", "tissue_ontology_term_id_2"],
                ),
            )

            response = self.app.post("/wmg/v1/query", json=request)

            self.assertEqual(200, response.status_code)

            # sanity check the expected values of the `n` stat for each data viz point; if this fails, the cube test
            # fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture has
            # changed
            expected_cell_count_per_cell_type = dim_size ** len(
                set(expression_summary_non_indexed_dims).difference({"cell_type_ontology_term_id"})
            )
            assert expected_cell_count_per_cell_type == 729

            # there are 729 possible combinations per tissue-cell type given the above filtering criteria,
            # and 10 cells per type in the cell counts cube, so we expect 7290 total cells per tissue-cell type

            expected_combinations_per_cell_type = dim_size ** len(
                set(expression_summary_non_indexed_dims).difference({"cell_type_ontology_term_id"})
            )
            expected_n_cells_per_cell_type = expected_combinations_per_cell_type * 10
            assert expected_n_cells_per_cell_type == 7290

            expected = {
                "snapshot_id": "dummy-snapshot",
                "expression_summary": {
                    "gene_ontology_term_id_0": {
                        "tissue_ontology_term_id_1": [
                            {
                                "id": "cell_type_ontology_term_id_0",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_1",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_2",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                        ],
                        "tissue_ontology_term_id_2": [
                            {
                                "id": "cell_type_ontology_term_id_0",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_1",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_2",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                        ],
                    },
                    "gene_ontology_term_id_2": {
                        "tissue_ontology_term_id_1": [
                            {
                                "id": "cell_type_ontology_term_id_0",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_1",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_2",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                        ],
                        "tissue_ontology_term_id_2": [
                            {
                                "id": "cell_type_ontology_term_id_0",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_1",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_2",
                                "n": 2187,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 2187 / (10 * (3**8)),
                            },
                        ],
                    },
                },
                "term_id_labels": {
                    "cell_types": {
                        "tissue_ontology_term_id_1": [
                            {
                                "cell_type": "cell_type_ontology_term_id_0_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                                "total_count": 21870,
                                "depth": 0,
                            },
                            {
                                "cell_type": "cell_type_ontology_term_id_1_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                                "total_count": 21870,
                                "depth": 1,
                            },
                            {
                                "cell_type": "cell_type_ontology_term_id_2_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                                "total_count": 21870,
                                "depth": 2,
                            },
                        ],
                        "tissue_ontology_term_id_2": [
                            {
                                "cell_type": "cell_type_ontology_term_id_0_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                                "total_count": 21870,
                                "depth": 0,
                            },
                            {
                                "cell_type": "cell_type_ontology_term_id_1_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                                "total_count": 21870,
                                "depth": 1,
                            },
                            {
                                "cell_type": "cell_type_ontology_term_id_2_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                                "total_count": 21870,
                                "depth": 2,
                            },
                        ],
                    },
                    "genes": [
                        {"gene_ontology_term_id_0": "gene_ontology_term_id_0_label"},
                        {"gene_ontology_term_id_2": "gene_ontology_term_id_2_label"},
                    ],
                },
            }
            self.assertEqual(expected, json.loads(response.data))

    @patch("backend.wmg.api.v1.gene_term_label")
    @patch("backend.wmg.api.v1.ontology_term_label")
    @patch("backend.wmg.api.v1.load_snapshot")
    def test__query_explicit_cell_ordering__returns_correct_cell_ordering(
        self, load_snapshot, ontology_term_label, gene_term_label
    ):
        dim_size = 2
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
            cell_ordering_generator_fn=reverse_cell_type_ordering,
        ) as snapshot:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_snapshot.return_value = snapshot

            # mock the functions in the ontology_labels module, so we can assert deterministic values in the
            # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
            # module is separately unit tested, and here we just want to verify the response building logic is correct.
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"

            request = dict(
                filter=dict(
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0", "tissue_ontology_term_id_1"],
                ),
            )

            response = self.app.post("/wmg/v1/query", json=request)

            self.assertEqual(200, response.status_code)

            expected = {
                "tissue_ontology_term_id_0": [
                    {
                        "cell_type": "cell_type_ontology_term_id_0_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                        "total_count": 1280,
                        "depth": 0,
                    },
                    {
                        "cell_type": "cell_type_ontology_term_id_1_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                        "total_count": 1280,
                        "depth": 1,
                    },
                ],
                "tissue_ontology_term_id_1": [
                    {
                        "cell_type": "cell_type_ontology_term_id_0_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                        "total_count": 1280,
                        "depth": 0,
                    },
                    {
                        "cell_type": "cell_type_ontology_term_id_1_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                        "total_count": 1280,
                        "depth": 1,
                    },
                ],
            }
            self.assertEqual(expected, json.loads(response.data)["term_id_labels"]["cell_types"])

    @patch("backend.wmg.api.v1.gene_term_label")
    @patch("backend.wmg.api.v1.ontology_term_label")
    @patch("backend.wmg.api.v1.load_snapshot")
    def test__query_total_cell_count_per_cell_type(self, load_snapshot, ontology_term_label, gene_term_label):
        expected_count = 42
        dim_size = 2
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=lambda coords: all_X_cell_counts_values(coords, expected_count),
            cell_ordering_generator_fn=reverse_cell_type_ordering,
        ) as snapshot:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_snapshot.return_value = snapshot

            # mock the functions in the ontology_labels module, so we can assert deterministic values in the
            # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
            # module is separately unit tested, and here we just want to verify the response building logic is correct.
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"

            request = dict(
                filter=dict(
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0", "tissue_ontology_term_id_1"],
                ),
            )

            response = self.app.post("/wmg/v1/query", json=request)

            self.assertEqual(200, response.status_code)

            # each cell type has `expected_count` cells for each possible combination of secondary filters
            # given the present constraints (1 organism, both tissues). There are 2**9=512 possible combinations
            # of filters. After aggregating the counts across two tissues and two cell types per tissue,
            # there are 128 entries per cell type-tissue combination. Hence, the toal count will be
            # expected_count * 128
            expected = {
                "tissue_ontology_term_id_0": [
                    {
                        "cell_type": "cell_type_ontology_term_id_0_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                        "total_count": expected_count * 128,
                        "depth": 0,
                    },
                    {
                        "cell_type": "cell_type_ontology_term_id_1_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                        "total_count": expected_count * 128,
                        "depth": 1,
                    },
                ],
                "tissue_ontology_term_id_1": [
                    {
                        "cell_type": "cell_type_ontology_term_id_0_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                        "total_count": expected_count * 128,
                        "depth": 0,
                    },
                    {
                        "cell_type": "cell_type_ontology_term_id_1_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                        "total_count": expected_count * 128,
                        "depth": 1,
                    },
                ],
            }
            self.assertEqual(expected, json.loads(response.data)["term_id_labels"]["cell_types"])

    def test__query_empty_request__returns_400(self):
        response = self.app.post("/wmg/v1/query", json={})

        self.assertEqual(400, response.status_code)
        self.assertEqual("'filter' is a required property", json.loads(response.data)["detail"])

    def test__query_missing_gene__request_returns_400(self):
        request = dict(
            filter=dict(
                organism_ontology_term_id="organism_ontology_term_id_0",
                tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            ),
        )

        response = self.app.post("/wmg/v1/query", json=request)

        self.assertEqual(400, response.status_code)

    def test__query_missing_organism__request_returns_400(self):
        request = dict(
            filter=dict(
                gene_ontology_term_ids=["gene_ontology_term_id_0"],
                tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            ),
        )

        response = self.app.post("/wmg/v1/query", json=request)

        self.assertEqual(400, response.status_code)

    def test__query_missing_tissue_request__returns_400(self):
        request = dict(
            filter=dict(
                gene_ontology_term_ids=["gene_ontology_term_id_0"],
                organism_ontology_term_id="organism_ontology_term_id_0",
            ),
        )

        response = self.app.post("/wmg/v1/query", json=request)

        self.assertEqual(400, response.status_code)

    """
    # these tests COULD work but the virtual test fixture is way too big to precompute filter relationships
    # because we cannot fit the cell counts dataframe in memory. These tests will be re-enabled once we
    # refine the test snapshot generation codepath.
    @patch("backend.wmg.api.v1.fetch_datasets_metadata")
    @patch("backend.wmg.api.v1.gene_term_label")
    @patch("backend.wmg.api.v1.ontology_term_label")
    @patch("backend.wmg.api.v1.load_snapshot")
    def test__query_request_with_filter_dims__returns_valid_filter_dims__base_case(
        self, load_snapshot, ontology_term_label, gene_term_label, fetch_datasets_metadata
    ):
        # mock the functions in the ontology_labels module, so we can assert deterministic values in the
        # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
        # module is separately unit tested, and here we just want to verify the response building logic is correct.
        dim_size = 1
        with create_temp_wmg_snapshot(dim_size=dim_size) as snapshot:
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"
            fetch_datasets_metadata.return_value = mock_datasets_metadata([f"dataset_id_{i}" for i in range(dim_size)])
            load_snapshot.return_value = snapshot
            filter_0 = dict(
                # these don't matter for the expected result
                gene_ontology_term_ids=["gene_ontology_term_id_0"],
                organism_ontology_term_id="organism_ontology_term_id_0",
                tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                dataset_ids=["dataset_id_0"],
                disease_ontology_term_ids=["disease_ontology_term_id_0"],
                sex_ontology_term_ids=["sex_ontology_term_id_0"],
                # these matter for the expected result
                development_stage_ontology_term_ids=["development_stage_ontology_term_id_0"],
                self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_0"],
            )

            filter_0_request = dict(
                filter=filter_0,
                # matters for this test
                is_rollup=True,
            )

            response = self.app.post("/wmg/v1/query", json=filter_0_request)
            expected_filters = {
                "datasets": [
                    {
                        "collection_id": "dataset_id_0_coll_id",
                        "collection_label": "dataset_id_0_coll_name",
                        "id": "dataset_id_0",
                        "label": "dataset_id_0_name",
                    }
                ],
                "development_stage_terms": [
                    {"development_stage_ontology_term_id_0": "development_stage_ontology_term_id_0_label"}
                ],
                "disease_terms": [{"disease_ontology_term_id_0": "disease_ontology_term_id_0_label"}],
                "self_reported_ethnicity_terms": [
                    {"self_reported_ethnicity_ontology_term_id_0": "self_reported_ethnicity_ontology_term_id_0_label"}
                ],
                "sex_terms": [{"sex_ontology_term_id_0": "sex_ontology_term_id_0_label"}],
                "tissue_terms": [{"tissue_ontology_term_id_0": "tissue_ontology_term_id_0_label"}],
            }
            self.assertEqual(json.loads(response.data)["filter_dims"], expected_filters)

    @patch("backend.wmg.api.v1.fetch_datasets_metadata")
    @patch("backend.wmg.api.v1.gene_term_label")
    @patch("backend.wmg.api.v1.ontology_term_label")
    @patch("backend.wmg.api.v1.load_snapshot")
    def test__query_request_with_filter_dims__returns_valid_filter_dims(
        self, load_snapshot, ontology_term_label, gene_term_label, fetch_datasets_metadata
    ):
        # mock the functions in the ontology_labels module, so we can assert deterministic values in the
        # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
        # module is separately unit tested, and here we just want to verify the response building logic is correct.
        dim_size = 3
        self.maxDiff = None
        with create_temp_wmg_snapshot(
            dim_size=dim_size, exclude_logical_coord_fn=exclude_dev_stage_and_ethnicity_for_secondary_filter_test
        ) as snapshot:
            # set up the expression summary cube for secondary filtering
            # drop all rows where ethnicity_1 and ethnicity_2 are associated with dev_stage_1 and dev_stage_2
            # thus filtering for dev_stage_0 should return filter options that include ethnicity 0,1 &2 but
            # filtering for dev_stage_1 or dev_stage_2 should only return ethnicity 0 (and vice versa)
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"
            fetch_datasets_metadata.return_value = mock_datasets_metadata([f"dataset_id_{i}" for i in range(dim_size)])
            # setup up API endpoints to use a mocked cube
            load_snapshot.return_value = snapshot
            with self.subTest("when a secondary dimension has criteria, it's own values are not restricted"):
                filter_0 = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    dataset_ids=["dataset_id_0"],
                    disease_ontology_term_ids=["disease_ontology_term_id_0"],
                    sex_ontology_term_ids=["sex_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=["development_stage_ontology_term_id_0"],
                    self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_0"],
                )

                filter_0_request = dict(
                    filter=filter_0,
                    # matters for this test
                    is_rollup=True,
                )

                filter_0_no_dev_stage_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    dataset_ids=["dataset_id_0"],
                    disease_ontology_term_ids=["disease_ontology_term_id_0"],
                    sex_ontology_term_ids=["sex_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=[],
                    self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_0"],
                )
                filter_0_no_dev_stage_request = dict(
                    filter=filter_0_no_dev_stage_filter, is_rollup=True
                )

                filter_0_no_ethnicity_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    dataset_ids=["dataset_id_0"],
                    disease_ontology_term_ids=["disease_ontology_term_id_0"],
                    sex_ontology_term_ids=["sex_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=["development_stage_ontology_term_id_0"],
                    self_reported_ethnicity_ontology_term_ids=[],
                )
                filter_0_no_ethnicity_request = dict(
                    filter=filter_0_no_ethnicity_filter, is_rollup=True
                )
                # the values for dev_stage terms when a dev stage filter is included should match the values returned
                # if no filter is passed in for dev stage
                response = self.app.post("/wmg/v1/query", json=filter_0_request)
                dev_stage_terms = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                self_reported_ethnicity_terms = json.loads(response.data)["filter_dims"][
                    "self_reported_ethnicity_terms"
                ]

                no_dev_stage_filter_response = self.app.post("/wmg/v1/query", json=filter_0_no_dev_stage_request)
                dev_stage_terms_if_no_dev_stage_filters = json.loads(no_dev_stage_filter_response.data)["filter_dims"][
                    "development_stage_terms"
                ]

                no_ethnicity_filter_response = self.app.post("/wmg/v1/query", json=filter_0_no_ethnicity_request)
                self_reported_ethnicity_terms_if_no_dev_stage_filters = json.loads(no_ethnicity_filter_response.data)[
                    "filter_dims"
                ]["self_reported_ethnicity_terms"]

                # filter options for dev_stage
                self.assertEqual(dev_stage_terms, dev_stage_terms_if_no_dev_stage_filters)
                self.assertEqual(self_reported_ethnicity_terms, self_reported_ethnicity_terms_if_no_dev_stage_filters)

            with self.subTest(
                "when a secondary dimension has criteria, the remaining filter values for other "
                "secondary dimensions are properly restricted"
            ):
                # filtering for dev_stage_0 should return all possible ethnicity terms
                # filtering for dev_stage_1 should only return ethnicity_0
                # filtering for dev_stage_2 should only return ethnicity_0
                # filtering for ethnicity_1 should only return dev_stage_0
                # filtering for ethnicity_2 should only return dev_stage_0

                # filtering for ethnicity_0 should return all dev stages
                # not filtering for dev should return all ethnicities
                self_reported_ethnicity_0_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=[],
                    self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_0"],
                )
                all_development_stage_terms = [
                    {"development_stage_ontology_term_id_0": "development_stage_ontology_term_id_0_label"},
                    {"development_stage_ontology_term_id_1": "development_stage_ontology_term_id_1_label"},
                    {"development_stage_ontology_term_id_2": "development_stage_ontology_term_id_2_label"},
                ]
                all_self_reported_ethnicity_terms = [
                    {"self_reported_ethnicity_ontology_term_id_0": "self_reported_ethnicity_ontology_term_id_0_label"},
                    {"self_reported_ethnicity_ontology_term_id_1": "self_reported_ethnicity_ontology_term_id_1_label"},
                    {"self_reported_ethnicity_ontology_term_id_2": "self_reported_ethnicity_ontology_term_id_2_label"},
                ]
                self_reported_ethnicity_0_request = dict(
                    filter=self_reported_ethnicity_0_filter, is_rollup=True
                )
                response = self.app.post("/wmg/v1/query", json=self_reported_ethnicity_0_request)
                dev_stage_terms = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                self_reported_ethnicity_terms = json.loads(response.data)["filter_dims"][
                    "self_reported_ethnicity_terms"
                ]
                self.assertEqual(all_development_stage_terms, dev_stage_terms)
                self.assertEqual(all_self_reported_ethnicity_terms, self_reported_ethnicity_terms)

                # filtering for ethnicity_1 should return dev_stage_0 (even when filtering for dev_stage_1 or 2)
                self_reported_ethnicity_1_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=[],
                    self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_1"],
                )
                # since secondary filters should not affect term options for their own category, filtering (or not
                # filtering) for dev stage 1 should not affect the terms returned
                expected_development_stage_terms = [
                    {"development_stage_ontology_term_id_0": "development_stage_ontology_term_id_0_label"},
                ]
                expected_self_reported_ethnicity_term = [
                    {"self_reported_ethnicity_ontology_term_id_0": "self_reported_ethnicity_ontology_term_id_0_label"}
                ]
                self_reported_ethnicity_1_request = dict(
                    filter=self_reported_ethnicity_1_filter, is_rollup=True
                )
                response = self.app.post("/wmg/v1/query", json=self_reported_ethnicity_1_request)
                dev_stage_terms_no_dev_filter = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                self_reported_ethnicity_terms_no_dev_filter = json.loads(response.data)["filter_dims"][
                    "self_reported_ethnicity_terms"
                ]
                self.assertEqual(expected_development_stage_terms, dev_stage_terms_no_dev_filter)
                self.assertEqual(all_self_reported_ethnicity_terms, self_reported_ethnicity_terms_no_dev_filter)

                self_reported_ethnicity_1_dev_1_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=["development_stage_ontology_term_id_1"],
                    self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_1"],
                )
                self_reported_ethnicity_1_dev_1_request = dict(
                    filter=self_reported_ethnicity_1_dev_1_filter, is_rollup=True
                )

                response = self.app.post("/wmg/v1/query", json=self_reported_ethnicity_1_dev_1_request)
                dev_stage_terms_dev_filter = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                self_reported_ethnicity_terms_dev_filter = json.loads(response.data)["filter_dims"][
                    "self_reported_ethnicity_terms"
                ]

                self.assertEqual(expected_development_stage_terms, dev_stage_terms_dev_filter)
                self.assertEqual(dev_stage_terms_no_dev_filter, dev_stage_terms_dev_filter)
                self.assertEqual(expected_self_reported_ethnicity_term, self_reported_ethnicity_terms_dev_filter)
                self.assertNotEqual(
                    self_reported_ethnicity_terms_no_dev_filter, self_reported_ethnicity_terms_dev_filter
                )

                # filtering for ethnicity_2 should return dev_stage_0 (even when filtering for dev_stage_1 or 2)
                self_reported_ethnicity_2_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=[],
                    self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_2"],
                )
                self_reported_ethnicity_2_dev_2_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=["development_stage_ontology_term_id_2"],
                    self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_2"],
                )
                self_reported_ethnicity_2_request = dict(
                    filter=self_reported_ethnicity_2_filter, is_rollup=True
                )
                eth_2_dev_2_request = dict(
                    filter=self_reported_ethnicity_2_dev_2_filter, is_rollup=True
                )
                response = self.app.post("/wmg/v1/query", json=self_reported_ethnicity_2_request)
                dev_stage_terms_eth_2_no_dev_filter = json.loads(response.data)["filter_dims"][
                    "development_stage_terms"
                ]
                self_reported_ethnicity_terms_eth_2_no_dev_filter = json.loads(response.data)["filter_dims"][
                    "self_reported_ethnicity_terms"
                ]
                self.assertEqual(expected_development_stage_terms, dev_stage_terms_eth_2_no_dev_filter)
                self.assertEqual(all_self_reported_ethnicity_terms, self_reported_ethnicity_terms_eth_2_no_dev_filter)

                response = self.app.post("/wmg/v1/query", json=eth_2_dev_2_request)
                dev_stage_terms_eth_2_dev_2 = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                eth_stage_terms_eth_2_dev_2 = json.loads(response.data)["filter_dims"]["self_reported_ethnicity_terms"]

                self.assertEqual(expected_development_stage_terms, dev_stage_terms_eth_2_dev_2)
                self.assertEqual(expected_self_reported_ethnicity_term, eth_stage_terms_eth_2_dev_2)
                self.assertEqual(dev_stage_terms_eth_2_dev_2, dev_stage_terms_eth_2_no_dev_filter)
                self.assertNotEqual(eth_stage_terms_eth_2_dev_2, self_reported_ethnicity_terms_eth_2_no_dev_filter)"""

    @patch("backend.wmg.api.v1.gene_term_label")
    @patch("backend.wmg.api.v1.ontology_term_label")
    @patch("backend.wmg.api.v1.load_snapshot")
    def test__markers_returns_200_and_correct_response(self, load_snapshot, ontology_term_label, gene_term_label):
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_snapshot.return_value = snapshot

            # mock the functions in the ontology_labels module, so we can assert deterministic values in the
            # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
            # module is separately unit tested, and here we just want to verify the response building logic is correct.
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"

            request = dict(
                celltype="CL:0000786",
                organism="NCBITaxon:9606",
                tissue="UBERON:0002048",
                n_markers=10,
                test="ttest",
            )

            response = self.app.post("/wmg/v1/markers", json=request)
            received = json.loads(response.data)

            expected = {
                "marker_genes": [
                    {"gene_ontology_term_id": "ENSG00000132465", "p_value": 0.0, "effect_size": 2.4539473056793213},
                    {"gene_ontology_term_id": "ENSG00000170476", "p_value": 0.0, "effect_size": 2.080190420150757},
                    {"gene_ontology_term_id": "ENSG00000180879", "p_value": 0.0, "effect_size": 2.0378074645996094},
                    {"gene_ontology_term_id": "ENSG00000134285", "p_value": 0.0, "effect_size": 1.7846676111221313},
                    {"gene_ontology_term_id": "ENSG00000099958", "p_value": 0.0, "effect_size": 1.574164628982544},
                    {"gene_ontology_term_id": "ENSG00000211592", "p_value": 0.0, "effect_size": 1.147048830986023},
                    {"gene_ontology_term_id": "ENSG00000166562", "p_value": 0.0, "effect_size": 1.1273812055587769},
                    {"gene_ontology_term_id": "ENSG00000118363", "p_value": 0.0, "effect_size": 1.0306953191757202},
                    {"gene_ontology_term_id": "ENSG00000125844", "p_value": 0.0, "effect_size": 0.46169212460517883},
                    {"gene_ontology_term_id": "ENSG00000100219", "p_value": 0.0, "effect_size": 0.4111901521682739},
                ],
                "snapshot_id": "realistic-test-snapshot",
            }
            self.assertDictEqual(received, expected)
            self.assertEqual(200, response.status_code)

    @patch("backend.wmg.api.v1.gene_term_label")
    @patch("backend.wmg.api.v1.ontology_term_label")
    @patch("backend.wmg.api.v1.load_snapshot")
    def test__markers_returns_200_and_empty_dictionary_for_bad_celltypes(
        self, load_snapshot, ontology_term_label, gene_term_label
    ):
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_snapshot.return_value = snapshot

            # mock the functions in the ontology_labels module, so we can assert deterministic values in the
            # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
            # module is separately unit tested, and here we just want to verify the response building logic is correct.
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"

            request = dict(
                celltype="CL:9999999",  # bad celltype
                organism="NCBITaxon:9606",
                tissue="UBERON:0002048",
                n_markers=10,
                test="ttest",
            )

            response = self.app.post("/wmg/v1/markers", json=request)
            received = json.loads(response.data)

            expected = {"marker_genes": [], "snapshot_id": "realistic-test-snapshot"}
            self.assertDictEqual(received, expected)
            self.assertEqual(200, response.status_code)


# mock the dataset and collection entity data that would otherwise be fetched from the db; in this test
# we only care that we're building the response correctly from the cube; WMG API integration tests verify
# with real datasets
def mock_datasets_metadata(dataset_ids):
    return [
        dict(
            id=dataset_id,
            label=f"{dataset_id}_name",
            collection_id=f"{dataset_id}_coll_id",
            collection_label=f"{dataset_id}_coll_name",
        )
        for dataset_id in dataset_ids
    ]
