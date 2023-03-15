import json
import unittest
from unittest.mock import patch

from backend.api_server.app import app
from backend.wmg.api.v1 import find_dimension_id_from_compare
from backend.wmg.data.query import MarkerGeneQueryCriteria
from tests.unit.backend.fixtures.environment_setup import EnvironmentSetup
from tests.unit.backend.wmg.fixtures.test_cube_schema import expression_summary_non_indexed_dims
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
    exclude_dev_stage_and_ethnicity_for_secondary_filter_test,
    load_realistic_test_snapshot,
    reverse_cell_type_ordering,
)
from tests.unit.backend.wmg.test_query import generate_expected_marker_gene_data_with_pandas

TEST_SNAPSHOT = "realistic-test-snapshot"


# this should only be used for generating expected outputs when using the test snapshot (see test_snapshot.py)
def generate_expected_term_id_labels_dictionary(genes, tissues, cell_types, total_count, compare_terms=None):
    result = {}
    result["cell_types"] = {}
    # assume tissues are sorted, and cell types are sorted within each tissue
    # assume the length of cell types is the dimensionality of each column in the
    # test snapshot (i.e. each tissue contains all fake cell types).
    # this line determines the starting order for each cell type in each tissue.
    orders = [int(tissue.split("_")[-1]) * len(cell_types) for tissue in tissues]
    for tissue, order in zip(tissues, orders):
        result["cell_types"][tissue] = {}
        for cell_type in cell_types:
            result["cell_types"][tissue][cell_type] = {}
            result["cell_types"][tissue][cell_type]["aggregated"] = {
                "cell_type_ontology_term_id": cell_type,
                "name": f"{cell_type}_label",
                "total_count": total_count,
                "order": order,
            }
            if compare_terms:
                for term in compare_terms:
                    result["cell_types"][tissue][cell_type][term] = {
                        "cell_type_ontology_term_id": cell_type,
                        "name": f"{term}_label",
                        "total_count": total_count // len(compare_terms),
                        "order": order,
                    }
            order += 1

    result["genes"] = []
    for gene in genes:
        result["genes"].append({gene: f"{gene}_label"})
    return result


def generate_expected_expression_summary_dictionary(genes, tissues, cell_types, n, me, pc, tpc, compare_terms=None):
    result = {}
    for gene in genes:
        result[gene] = {}
        for tissue in tissues:
            result[gene][tissue] = {}
            for cell_type in cell_types:
                result[gene][tissue][cell_type] = {}
                result[gene][tissue][cell_type]["aggregated"] = {
                    "n": n,
                    "me": me,
                    "pc": pc,
                    "tpc": tpc,
                }
                if compare_terms:
                    for term in compare_terms:
                        result[gene][tissue][cell_type][term] = {
                            "n": n // len(compare_terms),
                            "me": me,
                            "pc": pc,
                            "tpc": tpc / len(compare_terms),
                        }

    return result


def generate_test_inputs_and_expected_outputs(genes, tissues, organism, dim_size, me, expected_count, compare_dim=None):
    """
    Generates test inputs and expected outputs for the /wmg/v1/query endpoint.

    Arguments
    ---------
    genes: list of gene ontology term IDs
    tissues: list of tissue ontology term IDs
    organism: organism ontology term ID
    dim_size: size of each dimension of the test cube
    me: mean expression value to use for each gene/tissue/cell_type combination (scalar)
    expected_count: expected number of cells for each gene/tissue/cell_type combination (scalar)
    compare_dim: dimension to use for compare feature (optional). None if compare isn't used.

    Returns
    -------
    tuple of (request, expected_expression_summary, expected_term_id_labels)
    request: dictionary containing the request body to send to the /wmg/v1/query endpoint
    expected_expression_summary: dictionary containing the expected expression summary values
    expected_term_id_labels: dictionary containing the expected term ID labels
    """
    cell_types = [f"cell_type_ontology_term_id_{i}" for i in range(dim_size)]

    expected_combinations_per_cell_type = dim_size ** len(
        set(expression_summary_non_indexed_dims).difference({"cell_type_ontology_term_id"})
    )
    compare_terms = (
        [f"{find_dimension_id_from_compare(compare_dim)}_{i}" for i in range(dim_size)] if compare_dim else None
    )

    expected_term_id_labels = generate_expected_term_id_labels_dictionary(
        genes,
        tissues,
        cell_types,
        expected_combinations_per_cell_type * expected_count,
        compare_terms=compare_terms,
    )
    expected_expression_summary = generate_expected_expression_summary_dictionary(
        genes,
        tissues,
        cell_types,
        expected_combinations_per_cell_type,
        me,
        1 / expected_count,
        expected_combinations_per_cell_type / (expected_count * (dim_size ** len(expression_summary_non_indexed_dims))),
        compare_terms=compare_terms,
    )

    request = dict(
        filter=dict(
            gene_ontology_term_ids=genes,
            organism_ontology_term_id=organism,
            tissue_ontology_term_ids=tissues,
        )
    )
    if compare_dim:
        request["compare"] = compare_dim

    return (
        request,
        expected_expression_summary,
        expected_term_id_labels,
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

            genes = ["gene_ontology_term_id_0"]
            tissues = ["tissue_ontology_term_id_0"]
            organism = "organism_ontology_term_id_0"

            (request, expected_expression_summary, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, tissues, organism, dim_size, 1.0, 10
            )

            response = self.app.post("/wmg/v1/query", json=request)

            self.assertEqual(200, response.status_code)

            expected_response = {
                "snapshot_id": "dummy-snapshot",
                "expression_summary": expected_expression_summary,
                "term_id_labels": expected_term_id_labels,
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

            genes = ["gene_ontology_term_id_0", "gene_ontology_term_id_2"]
            tissues = ["tissue_ontology_term_id_1", "tissue_ontology_term_id_2"]
            organism = "organism_ontology_term_id_0"

            (request, expected_expression_summary, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, tissues, organism, dim_size, 1.0, 10
            )

            response = self.app.post("/wmg/v1/query", json=request)

            self.assertEqual(200, response.status_code)

            expected = {
                "snapshot_id": "dummy-snapshot",
                "expression_summary": expected_expression_summary,
                "term_id_labels": expected_term_id_labels,
            }
            self.assertEqual(expected, json.loads(response.data))

    @patch("backend.wmg.api.v1.gene_term_label")
    @patch("backend.wmg.api.v1.ontology_term_label")
    @patch("backend.wmg.api.v1.load_snapshot")
    def test__query_request_multi_primary_dims_only_with_compare__returns_200_and_correct_response(
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

            genes = ["gene_ontology_term_id_0", "gene_ontology_term_id_2"]
            tissues = ["tissue_ontology_term_id_1", "tissue_ontology_term_id_2"]
            organism = "organism_ontology_term_id_0"

            (request, expected_expression_summary, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, tissues, organism, dim_size, 1.0, 10, compare_dim="self_reported_ethnicity"
            )

            response = self.app.post("/wmg/v1/query", json=request)

            self.assertEqual(200, response.status_code)

            expected = {
                "snapshot_id": "dummy-snapshot",
                "expression_summary": expected_expression_summary,
                "term_id_labels": expected_term_id_labels,
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

            genes = ["gene_ontology_term_id_0"]
            tissues = ["tissue_ontology_term_id_0", "tissue_ontology_term_id_1"]
            organism = "organism_ontology_term_id_0"

            (request, _, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, tissues, organism, dim_size, 1.0, 10
            )

            response = self.app.post("/wmg/v1/query", json=request)

            self.assertEqual(200, response.status_code)

            expected = expected_term_id_labels["cell_types"]
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

            genes = ["gene_ontology_term_id_0"]
            tissues = ["tissue_ontology_term_id_0", "tissue_ontology_term_id_1"]
            organism = "organism_ontology_term_id_0"

            (request, _, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, tissues, organism, dim_size, 1.0, expected_count
            )

            response = self.app.post("/wmg/v1/query", json=request)

            self.assertEqual(200, response.status_code)

            expected = expected_term_id_labels["cell_types"]
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
                # these matter for the expected result
                development_stage_ontology_term_ids=["development_stage_ontology_term_id_0"],
                self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_0"],
            )

            filter_0_request = dict(
                filter=filter_0,
                is_rollup=True,
            )

            response = self.app.post("/wmg/v1/filters", json=filter_0_request)

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
                "disease_terms": [],
                "self_reported_ethnicity_terms": [
                    {"self_reported_ethnicity_ontology_term_id_0": "self_reported_ethnicity_ontology_term_id_0_label"}
                ],
                "sex_terms": [],
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
                    # these matter for the expected result
                    development_stage_ontology_term_ids=[],
                    self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_0"],
                )
                filter_0_no_dev_stage_request = dict(filter=filter_0_no_dev_stage_filter, is_rollup=True)

                filter_0_no_ethnicity_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=["development_stage_ontology_term_id_0"],
                    self_reported_ethnicity_ontology_term_ids=[],
                )
                filter_0_no_ethnicity_request = dict(filter=filter_0_no_ethnicity_filter, is_rollup=True)
                # the values for dev_stage terms when a dev stage filter is included should match the values returned
                # if no filter is passed in for dev stage
                response = self.app.post("/wmg/v1/filters", json=filter_0_request)
                dev_stage_terms = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                self_reported_ethnicity_terms = json.loads(response.data)["filter_dims"][
                    "self_reported_ethnicity_terms"
                ]

                no_dev_stage_filter_response = self.app.post("/wmg/v1/filters", json=filter_0_no_dev_stage_request)
                dev_stage_terms_if_no_dev_stage_filters = json.loads(no_dev_stage_filter_response.data)["filter_dims"][
                    "development_stage_terms"
                ]

                no_ethnicity_filter_response = self.app.post("/wmg/v1/filters", json=filter_0_no_ethnicity_request)
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
                self_reported_ethnicity_0_request = dict(filter=self_reported_ethnicity_0_filter, is_rollup=True)
                response = self.app.post("/wmg/v1/filters", json=self_reported_ethnicity_0_request)
                dev_stage_terms = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["development_stage_terms"]
                )
                self_reported_ethnicity_terms = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["self_reported_ethnicity_terms"]
                )
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
                self_reported_ethnicity_1_request = dict(filter=self_reported_ethnicity_1_filter, is_rollup=True)
                response = self.app.post("/wmg/v1/filters", json=self_reported_ethnicity_1_request)
                dev_stage_terms_no_dev_filter = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["development_stage_terms"]
                )
                self_reported_ethnicity_terms_no_dev_filter = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["self_reported_ethnicity_terms"]
                )
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

                response = self.app.post("/wmg/v1/filters", json=self_reported_ethnicity_1_dev_1_request)
                dev_stage_terms_dev_filter = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["development_stage_terms"]
                )
                self_reported_ethnicity_terms_dev_filter = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["self_reported_ethnicity_terms"]
                )

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
                self_reported_ethnicity_2_request = dict(filter=self_reported_ethnicity_2_filter, is_rollup=True)
                eth_2_dev_2_request = dict(filter=self_reported_ethnicity_2_dev_2_filter, is_rollup=True)
                response = self.app.post("/wmg/v1/filters", json=self_reported_ethnicity_2_request)
                dev_stage_terms_eth_2_no_dev_filter = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["development_stage_terms"]
                )
                self_reported_ethnicity_terms_eth_2_no_dev_filter = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["self_reported_ethnicity_terms"]
                )
                self.assertEqual(expected_development_stage_terms, dev_stage_terms_eth_2_no_dev_filter)
                self.assertEqual(all_self_reported_ethnicity_terms, self_reported_ethnicity_terms_eth_2_no_dev_filter)

                response = self.app.post("/wmg/v1/filters", json=eth_2_dev_2_request)
                dev_stage_terms_eth_2_dev_2 = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["development_stage_terms"]
                )
                eth_stage_terms_eth_2_dev_2 = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["self_reported_ethnicity_terms"]
                )

                self.assertEqual(expected_development_stage_terms, dev_stage_terms_eth_2_dev_2)
                self.assertEqual(expected_self_reported_ethnicity_term, eth_stage_terms_eth_2_dev_2)
                self.assertEqual(dev_stage_terms_eth_2_dev_2, dev_stage_terms_eth_2_no_dev_filter)
                self.assertNotEqual(eth_stage_terms_eth_2_dev_2, self_reported_ethnicity_terms_eth_2_no_dev_filter)

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
            criteria = MarkerGeneQueryCriteria(
                tissue_ontology_term_id=request["tissue"],
                cell_type_ontology_term_id=request["celltype"],
                organism_ontology_term_id=request["organism"],
            )
            expected = {
                "marker_genes": generate_expected_marker_gene_data_with_pandas(snapshot, criteria, "ttest", 10),
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


def _sort_list_by_dictionary_keys(d):
    return sorted(d, key=lambda k: list(k.keys())[0])
