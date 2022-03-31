import json
import unittest
from typing import Tuple
from unittest.mock import patch

from backend.corpora.api_server.app import app
from backend.wmg.data.schemas.cube_schema import cube_non_indexed_dims
from tests.unit.backend.corpora.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.wmg.fixtures.test_cube import (
    create_temp_wmg_cubes,
    all_ones_expression_summary_values,
    all_tens_cell_counts_values,
)


@unittest.skip("TileDB bug (<=0.13.1) causing these to fail")
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

    @patch("backend.wmg.api.v1.load_cubes")
    def test__primary_filter_dimensions__returns_200(self, load_cubes):
        # This test appears to be hitting a TileDB (<=0.13.1) bug and fails (intermittently) if dim_size=3
        with create_temp_wmg_cubes(dim_size=1) as cubes:
            # setup up API endpoints to use a mocked cube
            load_cubes.return_value = cubes

            response = self.app.get("/wmg/v1/primary_filter_dimensions")

        self.assertEqual(200, response.status_code)

    @patch("backend.wmg.data.query.gene_term_label")
    @patch("backend.wmg.data.query.ontology_term_label")
    @patch("backend.wmg.api.v1.load_cubes")
    def test__primary_filter_dimensions__returns_valid_response_body(
        self, load_cubes, ontology_term_label, gene_term_label
    ):

        # we want disjoint set of genes across organisms, to mimic reality (each organism has its own set of genes);
        # without this filtering function, the cube would have the cross-product of organisms * genes
        def exclude(logical_coord: Tuple) -> bool:
            return (logical_coord[0], logical_coord[2]) not in {
                ("gene_ontology_term_id_0", "organism_ontology_term_id_0"),
                ("gene_ontology_term_id_1", "organism_ontology_term_id_1"),
                ("gene_ontology_term_id_2", "organism_ontology_term_id_2"),
            }

        # This test appears to be hitting a TileDB (<=0.13.1) bug and fails if dim_size=3
        with create_temp_wmg_cubes(exclude_logical_coord_fn=exclude, dim_size=1) as cubes:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_cubes.return_value = cubes

            # mock the functions in the ontology_labels module, so we can assert deterministic values in the
            # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
            # module is separately unit tested, and here we just want to verify the response building logic is correct.
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"

            response = self.app.get("/wmg/v1/primary_filter_dimensions")

        expected = dict(
            snapshot_id="dummy-snapshot",
            organism_terms=[
                {"organism_ontology_term_id_0": "organism_ontology_term_id_0_label"},
                # {"organism_ontology_term_id_1": "organism_ontology_term_id_1_label"},
                # {"organism_ontology_term_id_2": "organism_ontology_term_id_2_label"},
            ],
            tissue_terms=[
                {"tissue_ontology_term_id_0": "tissue_ontology_term_id_0_label"},
                # {"tissue_ontology_term_id_1": "tissue_ontology_term_id_1_label"},
                # {"tissue_ontology_term_id_2": "tissue_ontology_term_id_2_label"},
            ],
            gene_terms={
                "organism_ontology_term_id_0": [
                    {"gene_ontology_term_id_0": "gene_ontology_term_id_0_label"},
                ],
                # "organism_ontology_term_id_1": [
                #     {"gene_ontology_term_id_1": "gene_ontology_term_id_1_label"},
                # ],
                # "organism_ontology_term_id_2": [{"gene_ontology_term_id_2": "gene_ontology_term_id_2_label"}],
            },
        )

        self.assertEqual(expected, json.loads(response.data))

    @patch("backend.wmg.data.query.gene_term_label")
    @patch("backend.wmg.data.query.ontology_term_label")
    @patch("backend.wmg.api.v1.load_cubes")
    def test__query_minimal_valid_request__returns_200_and_empty_expr_summary(
        self, load_cubes, ontology_term_label, gene_term_label
    ):
        dim_size = 1
        with create_temp_wmg_cubes(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as cubes:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_cubes.return_value = cubes

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
                    "cell_types": [{"cell_type_ontology_term_id_0": "cell_type_ontology_term_id_0_label"}],
                    "genes": [{"gene_ontology_term_id_0": "gene_ontology_term_id_0_label"}],
                },
                "filter_dims": {},
            }
            self.assertEqual(expected_response, json.loads(response.data))

    def test__query_empty_request__returns_400(self):
        response = self.app.post("/wmg/v1/query", json={})

        self.assertEqual(400, response.status_code)
        self.assertEqual("'filter' is a required property", json.loads(response.data)["detail"])

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

    @patch("backend.wmg.data.query.gene_term_label")
    @patch("backend.wmg.data.query.ontology_term_label")
    @patch("backend.wmg.api.v1.load_cubes")
    def test__query_request_primary_dims_only__returns_valid_response_body(
        self, load_cubes, ontology_term_label, gene_term_label
    ):
        dim_size = 3
        with create_temp_wmg_cubes(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as cubes:
            # setup up API endpoints to use a mocked cube containing all stat values of 1, for a deterministic
            # expected query response
            load_cubes.return_value = cubes

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
                "filter_dims": {},
            }
            self.assertEqual(expected, json.loads(response.data))

    @patch("backend.wmg.api.v1.fetch_datasets_metadata")
    @patch("backend.wmg.data.query.gene_term_label")
    @patch("backend.wmg.data.query.ontology_term_label")
    @patch("backend.wmg.api.v1.load_cubes")
    def test__query_request_with_filter_dims__returns_valid_filter_dims(
        self, load_cubes, ontology_term_label, gene_term_label, fetch_datasets_metadata
    ):
        dim_size = 3
        with create_temp_wmg_cubes(dim_size=dim_size) as cubes:
            # setup up API endpoints to use a mocked cube
            load_cubes.return_value = cubes

            # mock the functions in the ontology_labels module, so we can assert deterministic values in the
            # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
            # module is separately unit tested, and here we just want to verify the response building logic is correct.
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"

            fetch_datasets_metadata.return_value = mock_datasets_metadata([f"dataset_id_{i}" for i in range(dim_size)])

            request = dict(
                # doesn't matter for this test
                filter=dict(
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                ),
                # matters for this test
                include_filter_dims=True,
            )

            response = self.app.post("/wmg/v1/query", json=request)

            expected = {
                "datasets": [
                    {
                        "id": "dataset_id_0",
                        "label": "dataset_id_0_name",
                        "collection_label": "dataset_id_0_coll_name",
                        "collection_id": "dataset_id_0_coll_id",
                    },
                    {
                        "id": "dataset_id_1",
                        "label": "dataset_id_1_name",
                        "collection_label": "dataset_id_1_coll_name",
                        "collection_id": "dataset_id_1_coll_id",
                    },
                    {
                        "id": "dataset_id_2",
                        "label": "dataset_id_2_name",
                        "collection_label": "dataset_id_2_coll_name",
                        "collection_id": "dataset_id_2_coll_id",
                    },
                ],
                "development_stage_terms": [
                    {"development_stage_ontology_term_id_0": "development_stage_ontology_term_id_0_label"},
                    {"development_stage_ontology_term_id_1": "development_stage_ontology_term_id_1_label"},
                    {"development_stage_ontology_term_id_2": "development_stage_ontology_term_id_2_label"},
                ],
                "disease_terms": [
                    {"disease_ontology_term_id_0": "disease_ontology_term_id_0_label"},
                    {"disease_ontology_term_id_1": "disease_ontology_term_id_1_label"},
                    {"disease_ontology_term_id_2": "disease_ontology_term_id_2_label"},
                ],
                "ethnicity_terms": [
                    {"ethnicity_ontology_term_id_0": "ethnicity_ontology_term_id_0_label"},
                    {"ethnicity_ontology_term_id_1": "ethnicity_ontology_term_id_1_label"},
                    {"ethnicity_ontology_term_id_2": "ethnicity_ontology_term_id_2_label"},
                ],
                "sex_terms": [
                    {"sex_ontology_term_id_0": "sex_ontology_term_id_0_label"},
                    {"sex_ontology_term_id_1": "sex_ontology_term_id_1_label"},
                    {"sex_ontology_term_id_2": "sex_ontology_term_id_2_label"},
                ],
            }
            self.assertEqual(expected, json.loads(response.data)["filter_dims"])

    @patch("backend.wmg.api.v1.fetch_datasets_metadata")
    @patch("backend.wmg.data.query.gene_term_label")
    @patch("backend.wmg.data.query.ontology_term_label")
    @patch("backend.wmg.api.v1.load_cubes")
    def test__query_request_with_filter_dims__returns_elided_filter_dims(
        self, load_cubes, ontology_term_label, gene_term_label, fetch_datasets_metadata
    ):
        dim_size = 2
        with create_temp_wmg_cubes(dim_size=dim_size) as cubes:
            load_cubes.return_value = cubes

            # mock the functions in the ontology_labels module, so we can assert deterministic values in the
            # "term_id_labels" portion of the response body; note that the correct behavior of the ontology_labels
            # module is separately unit tested, and here we just want to verify the response building logic is correct.
            ontology_term_label.side_effect = lambda ontology_term_id: f"{ontology_term_id}_label"
            gene_term_label.side_effect = lambda gene_term_id: f"{gene_term_id}_label"

            fetch_datasets_metadata.return_value = mock_datasets_metadata(["dataset_id_0"])

            request = dict(
                filter=dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    dataset_ids=["dataset_id_0"],
                    disease_ontology_term_ids=["disease_ontology_term_id_1"],
                    sex_ontology_term_ids=["sex_ontology_term_id_0"],
                    development_stage_ontology_term_ids=["development_stage_ontology_term_id_1"],
                    ethnicity_ontology_term_ids=["ethnicity_ontology_term_id_0"],
                ),
                include_filter_dims=True,
            )

            response = self.app.post("/wmg/v1/query", json=request)

            expected_filter_dims = {
                "datasets": [
                    {
                        "id": "dataset_id_0",
                        "label": "dataset_id_0_name",
                        "collection_label": "dataset_id_0_coll_name",
                        "collection_id": "dataset_id_0_coll_id",
                    },
                ],
                "disease_terms": [
                    {"disease_ontology_term_id_1": "disease_ontology_term_id_1_label"},
                ],
                "sex_terms": [
                    {"sex_ontology_term_id_0": "sex_ontology_term_id_0_label"},
                ],
                "development_stage_terms": [
                    {"development_stage_ontology_term_id_1": "development_stage_ontology_term_id_1_label"},
                ],
                "ethnicity_terms": [
                    {"ethnicity_ontology_term_id_0": "ethnicity_ontology_term_id_0_label"},
                ],
            }
            self.maxDiff = None
            self.assertEqual(expected_filter_dims, json.loads(response.data)["filter_dims"])


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
