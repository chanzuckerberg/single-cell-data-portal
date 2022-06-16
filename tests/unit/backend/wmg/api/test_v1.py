import json
import unittest
from unittest.mock import patch

from backend.corpora.api_server.app import app
from backend.wmg.data.schemas.cube_schema import cube_non_indexed_dims
from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.wmg.fixtures.test_primary_filters import (
    test_snapshot_id,
    test_organism_terms,
    test_tissue_terms,
    test_gene_terms,
)
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    create_temp_wmg_snapshot,
    all_ones_expression_summary_values,
    all_tens_cell_counts_values,
    reverse_cell_type_ordering,
    exclude_all_but_one_gene_per_organism,
)


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
                set(cube_non_indexed_dims).difference({"cell_type_ontology_term_id"})
            )
            assert expected_cell_count_per_cell_type == 729

            expected = {
                "snapshot_id": "dummy-snapshot",
                "expression_summary": {
                    "gene_ontology_term_id_0": {
                        "tissue_ontology_term_id_1": [
                            {
                                "id": "cell_type_ontology_term_id_0",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_1",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_2",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                        ],
                        "tissue_ontology_term_id_2": [
                            {
                                "id": "cell_type_ontology_term_id_0",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_1",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_2",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                        ],
                    },
                    "gene_ontology_term_id_2": {
                        "tissue_ontology_term_id_1": [
                            {
                                "id": "cell_type_ontology_term_id_0",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_1",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_2",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                        ],
                        "tissue_ontology_term_id_2": [
                            {
                                "id": "cell_type_ontology_term_id_0",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_1",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
                            },
                            {
                                "id": "cell_type_ontology_term_id_2",
                                "n": 729,
                                "me": 1.0,
                                "pc": 0.1,
                                "tpc": 729 / (10 * (3**7)),
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
                                "depth": 0,
                            },
                            {
                                "cell_type": "cell_type_ontology_term_id_1_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                                "depth": 1,
                            },
                            {
                                "cell_type": "cell_type_ontology_term_id_2_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                                "depth": 2,
                            },
                        ],
                        "tissue_ontology_term_id_2": [
                            {
                                "cell_type": "cell_type_ontology_term_id_0_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                                "depth": 0,
                            },
                            {
                                "cell_type": "cell_type_ontology_term_id_1_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                                "depth": 1,
                            },
                            {
                                "cell_type": "cell_type_ontology_term_id_2_label",
                                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                                "depth": 2,
                            },
                        ],
                    },
                    "genes": [
                        {"gene_ontology_term_id_0": "gene_ontology_term_id_0_label"},
                        {"gene_ontology_term_id_2": "gene_ontology_term_id_2_label"},
                    ],
                },
                "filter_dims": {},
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
                        "depth": 0,
                    },
                    {
                        "cell_type": "cell_type_ontology_term_id_1_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                        "depth": 1,
                    },
                ],
                "tissue_ontology_term_id_1": [
                    {
                        "cell_type": "cell_type_ontology_term_id_0_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                        "depth": 0,
                    },
                    {
                        "cell_type": "cell_type_ontology_term_id_1_label",
                        "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
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
        ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
        gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"
        fetch_datasets_metadata.return_value = mock_datasets_metadata([f"dataset_id_{i}" for i in range(dim_size)])

        with create_temp_wmg_snapshot(dim_size=dim_size) as snapshot:
            # set up the expression summary cube for secondary filtering
            # drop all rows where ethnicity_1 and ethnicity_2 are associated with dev_stage_1 and dev_stage_2
            # thus filtering for dev_stage_0 should return filter options that include ethnicity 0,1 &2 but
            # filtering for dev_stage_1 or dev_stage_2 should only return ethnicity 0 (and vice versa)

            df = snapshot.expression_summary_cube.df[:]

            df.drop(
                df[
                    (
                        df.development_stage_ontology_term_id.isin(
                            ["development_stage_ontology_term_id_1", "development_stage_ontology_term_id_2"]
                        )
                    )
                    & (
                        df.ethnicity_ontology_term_id.isin(
                            ["ethnicity_ontology_term_id_1", "ethnicity_ontology_term_id_2"]
                        )
                    )
                ].index,
                inplace=True,
            )
            # check that rows were dropped
            self.assertGreater(snapshot.expression_summary_cube.df[:].shape[0], df.shape[0])
            # setup up API endpoints to use a mocked cube
            snapshot.expression_summary_cube.df[:] = df
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
                    ethnicity_ontology_term_ids=["ethnicity_ontology_term_id_0"],
                )

                filter_0_request = dict(
                    filter=filter_0,
                    # matters for this test
                    include_filter_dims=True,
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
                    ethnicity_ontology_term_ids=["ethnicity_ontology_term_id_0"],
                )
                filter_0_no_dev_stage_request = dict(filter=filter_0_no_dev_stage_filter, include_filter_dims=True)

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
                    ethnicity_ontology_term_ids=[],
                )
                filter_0_no_ethnicity_request = dict(filter=filter_0_no_ethnicity_filter, include_filter_dims=True)
                # the values for dev_stage terms when a dev stage filter is included should match the values returned
                # if no filter is passed in for dev stage
                response = self.app.post("/wmg/v1/query", json=filter_0_request)
                dev_stage_terms = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                ethnicity_terms = json.loads(response.data)["filter_dims"]["ethnicity_terms"]

                no_dev_stage_filter_response = self.app.post("/wmg/v1/query", json=filter_0_no_dev_stage_request)
                dev_stage_terms_if_no_dev_stage_filters = json.loads(no_dev_stage_filter_response.data)["filter_dims"][
                    "development_stage_terms"
                ]

                no_ethnicity_filter_response = self.app.post("/wmg/v1/query", json=filter_0_no_ethnicity_request)
                ethnictiy_terms_if_no_dev_stage_filters = json.loads(no_ethnicity_filter_response.data)["filter_dims"][
                    "development_stage_terms"
                ]

                # filter options for dev_stage
                self.assertEqual(dev_stage_terms, dev_stage_terms_if_no_dev_stage_filters)
                self.assertEqual(ethnicity_terms, ethnictiy_terms_if_no_dev_stage_filters)

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
                ethnicity_0_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=[],
                    ethnicity_ontology_term_ids=["ethnicity_ontology_term_id_0"],
                )
                all_development_stage_terms = [
                    {"development_stage_ontology_term_id_0": "development_stage_ontology_term_id_0_label"},
                    {"development_stage_ontology_term_id_1": "development_stage_ontology_term_id_1_label"},
                    {"development_stage_ontology_term_id_2": "development_stage_ontology_term_id_2_label"},
                ]
                all_ethnicity_terms = [
                    {"ethnicity_ontology_term_id_0": "ethnicity_ontology_term_id_0_label"},
                    {"ethnicity_ontology_term_id_1": "ethnicity_ontology_term_id_1_label"},
                    {"ethnicity_ontology_term_id_2": "ethnicity_ontology_term_id_2_label"},
                ]
                ethnicity_0_request = dict(filter=ethnicity_0_filter, include_filter_dims=True)
                response = self.app.post("/wmg/v1/query", json=ethnicity_0_request)
                dev_stage_terms = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                ethnicity_terms = json.loads(response.data)["filter_dims"]["ethnicity_terms"]
                self.assertEqual(all_development_stage_terms, dev_stage_terms)
                self.assertEqual(all_ethnicity_terms, ethnicity_terms)

                # filtering for ethnicity_1 should return dev_stage_0 (even when filtering for dev_stage_1 or 2)
                ethnicity_1_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=[],
                    ethnicity_ontology_term_ids=["ethnicity_ontology_term_id_1"],
                )
                # since secondary filters should not affect term options for their own category, filtering (or not
                # filtering) for dev stage 1 should not affect the terms returned
                expected_development_stage_terms = [
                    {"development_stage_ontology_term_id_0": "development_stage_ontology_term_id_0_label"},
                ]
                expected_ethnicity_term = all_ethnicity_terms = [
                    {"ethnicity_ontology_term_id_0": "ethnicity_ontology_term_id_0_label"}
                ]

                ethnicity_1_request = dict(filter=ethnicity_1_filter, include_filter_dims=True)
                response = self.app.post("/wmg/v1/query", json=ethnicity_1_request)
                dev_stage_terms_no_dev_filter = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                ethnicity_terms_no_dev_filter = json.loads(response.data)["filter_dims"]["ethnicity_terms"]
                self.assertEqual(expected_development_stage_terms, dev_stage_terms_no_dev_filter)
                self.assertEqual(all_ethnicity_terms, ethnicity_terms_no_dev_filter)

                ethnicity_1_dev_1_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=["development_stage_ontology_term_id_1"],
                    ethnicity_ontology_term_ids=["ethnicity_ontology_term_id_1"],
                )
                ethnicity_1_dev_1_request = dict(filter=ethnicity_1_dev_1_filter, include_filter_dims=True)
                response = self.app.post("/wmg/v1/query", json=ethnicity_1_dev_1_request)
                dev_stage_terms_dev_filter = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                ethnicity_terms_dev_filter = json.loads(response.data)["filter_dims"]["ethnicity_terms"]

                self.assertEqual(expected_development_stage_terms, dev_stage_terms_dev_filter)
                self.assertEqual(dev_stage_terms_no_dev_filter, dev_stage_terms_dev_filter)
                self.assertEqual(expected_ethnicity_term, ethnicity_terms_dev_filter)
                self.assertNotEqual(ethnicity_terms_no_dev_filter, ethnicity_terms_dev_filter)

                # filtering for ethnicity_2 should return dev_stage_0 (even when filtering for dev_stage_1 or 2)
                ethnicity_2_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=[],
                    ethnicity_ontology_term_ids=["ethnicity_ontology_term_id_2"],
                )
                ethnicity_2_dev_2_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=["development_stage_ontology_term_id_2"],
                    ethnicity_ontology_term_ids=["ethnicity_ontology_term_id_2"],
                )
                ethnicity_2_request = dict(filter=ethnicity_2_filter, include_filter_dims=True)
                eth_2_dev_2_request = dict(filter=ethnicity_2_dev_2_filter, include_filter_dims=True)
                response = self.app.post("/wmg/v1/query", json=ethnicity_2_request)
                dev_stage_terms_eth_2_no_dev_filter = json.loads(response.data)["filter_dims"][
                    "development_stage_terms"
                ]
                ethnicity_terms_eth_2_no_dev_filter = json.loads(response.data)["filter_dims"]["ethnicity_terms"]
                self.assertEqual(expected_development_stage_terms, dev_stage_terms_eth_2_no_dev_filter)
                self.assertEqual(all_ethnicity_terms, ethnicity_terms_eth_2_no_dev_filter)

                response = self.app.post("/wmg/v1/query", json=eth_2_dev_2_request)
                dev_stage_terms_eth_2_dev_2 = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                eth_stage_terms_eth_2_dev_2 = json.loads(response.data)["filter_dims"]["ethnicity_terms"]

                self.assertEqual(expected_development_stage_terms, dev_stage_terms_eth_2_dev_2)
                self.assertEqual(expected_ethnicity_term, eth_stage_terms_eth_2_dev_2)
                self.assertEqual(dev_stage_terms_eth_2_dev_2, dev_stage_terms_eth_2_no_dev_filter)
                self.assertNotEqual(eth_stage_terms_eth_2_dev_2, ethnicity_terms_eth_2_no_dev_filter)


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
