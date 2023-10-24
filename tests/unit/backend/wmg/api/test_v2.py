import json
import unittest
from unittest.mock import patch

from pytest import approx

from backend.api_server.app import app
from backend.wmg.api.v2 import find_dimension_id_from_compare
from backend.wmg.data.query import MarkerGeneQueryCriteria
from tests.test_utils import compare_dicts
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
    forward_cell_type_ordering,
    load_realistic_test_snapshot,
    reverse_cell_type_ordering,
)
from tests.unit.backend.wmg.test_query import generate_expected_marker_gene_data_with_pandas

TEST_SNAPSHOT = "realistic-test-snapshot"


# this should only be used for generating expected outputs when using the test snapshot (see test_snapshot.py)
def generate_expected_term_id_labels_dictionary(
    *,
    genes: list[str],
    tissues: list[str],
    cell_types: list[str],
    cell_count_tissue_cell_type: int,
    compare_terms: list[str],
    cell_counts_tissue_cell_type_compare_dim: int,
    cell_ordering_func=forward_cell_type_ordering,
) -> dict:
    """
    Generates aggregated cell counts and cell ordering expected to be returned by /wmg/v2/query endpoint.

    Arguments
    ---------
    genes: list of gene ontology term IDs

    tissues: list of tissue ontology term IDs

    cell_types: list of cell_type ontology term IDs

    cell_count_tissue_cell_type: total number of cells for each (tissue, cell_type) combination (scalar)

    compare_terms: list of ontology term IDs form the compare dimension (optional)

    cell_counts_tissue_cell_type_compare_dim: total number of cells for each
    (tissue, cell_type, <compare_dim>) combination (scalar)

    Returns
    -------
    result: A dictionary containing aggregate cell counts and cell type ordering info that is indexable by
            ["cell_types"][<tissue_ontology_term_id>][<cell_type_ontology_term_id>][<compare_dimension>]
    """

    result = {}
    result["cell_types"] = {}

    orders = sum([cell_ordering_func(cell_types) for _ in tissues], [])
    index = 0
    for tissue in tissues:
        result["cell_types"][tissue] = {}

        for cell_type in cell_types:
            result["cell_types"][tissue][cell_type] = {}
            result["cell_types"][tissue][cell_type]["aggregated"] = {
                "cell_type_ontology_term_id": cell_type,
                "name": f"{cell_type}_label",
                "total_count": cell_count_tissue_cell_type,
                "order": orders[index],
            }

            for term in compare_terms:
                result["cell_types"][tissue][cell_type][term] = {
                    "cell_type_ontology_term_id": cell_type,
                    "name": f"{term}_label",
                    "total_count": cell_counts_tissue_cell_type_compare_dim,
                    "order": orders[index],
                }
            index += 1

        tissue_cell_counts = {
            "tissue_ontology_term_id": tissue,
            "name": f"{tissue}_label",
            "total_count": sum(
                [agg_dict["aggregated"]["total_count"] for _, agg_dict in result["cell_types"][tissue].items()]
            ),
            "order": order,
        }

        result["cell_types"][tissue]["tissue_stats"] = {}
        result["cell_types"][tissue]["tissue_stats"]["aggregated"] = tissue_cell_counts

    result["genes"] = []
    for gene in genes:
        result["genes"].append({gene: f"{gene}_label"})

    return result


def generate_expected_expression_summary_dictionary(
    *,
    genes: list[str],
    tissues: list[str],
    cell_count_tissue: int,
    cell_types: list[str],
    cell_count_tissue_cell_type: int,
    nnz_gene_tissue_cell_type: int,
    compare_terms: list[str],
    cell_counts_tissue_cell_type_compare_dim: int,
    nnz_gene_tissue_cell_type_compare_dim: int,
    me: float,
) -> dict:
    """
    Generates expression summary stats expected to be returned by /wmg/v2/query endpoint.

    Arguments
    ---------
    genes: list of gene ontology term IDs

    tissues: list of tissue ontology term IDs

    cell_count_tissue: total number of cells for each tissue combination (scalar)

    cell_types: list of cell_type ontology term IDs

    cell_count_tissue_cell_type: total number of cells for each (tissue, cell_type) combination (scalar)

    nnz_gene_tissue_cell_type: the nnz value for each (gene, tissue, cell_type) combination (scalar)

    compare_terms: list of ontology term IDs form the compare dimension (optional)

    cell_counts_tissue_cell_type_compare_dim: total number of cells for each
    (tissue, cell_type, <compare_dim>) combination (scalar)

    nnz_gene_tissue_cell_type_compare_dim: the nnz value for each (gene, tissue, cell_type, <compare_dim>) combination (scalar)

    me: mean expression value to use for each (gene, tissue, cell_type) combination (scalar)

    Returns
    -------
    result: A dictionary containing gene expression stats that is indexable by
            [<gene_ontology_term_id>][<tissue_ontology_term_id>][<cell_type_ontology_term_id>][<compare_dimension>]
    """
    result = {}
    for gene in genes:
        if gene != ".":
            result[gene] = {}
            for tissue in tissues:
                result[gene][tissue] = {}

                for cell_type in cell_types:
                    result[gene][tissue][cell_type] = {}

                    pc_gene_tissue_cell_type = nnz_gene_tissue_cell_type / cell_count_tissue_cell_type
                    tpc_gene_tissue_cell_type = nnz_gene_tissue_cell_type / cell_count_tissue

                    result[gene][tissue][cell_type]["aggregated"] = {
                        "n": nnz_gene_tissue_cell_type,
                        "me": me,
                        "pc": pc_gene_tissue_cell_type,
                        "tpc": tpc_gene_tissue_cell_type,
                    }

                    for term in compare_terms:
                        pc_gene_tissue_cell_type_compare_dim = (
                            nnz_gene_tissue_cell_type_compare_dim / cell_counts_tissue_cell_type_compare_dim
                        )
                        tpc_gene_tissue_cell_type_compare_dim = (
                            nnz_gene_tissue_cell_type_compare_dim / cell_count_tissue
                        )

                        result[gene][tissue][cell_type][term] = {
                            "n": nnz_gene_tissue_cell_type_compare_dim,
                            "me": me,
                            "pc": pc_gene_tissue_cell_type_compare_dim,
                            "tpc": tpc_gene_tissue_cell_type_compare_dim,
                        }

                nnz_gene_tissue = sum([agg_dict["aggregated"]["n"] for _, agg_dict in result[gene][tissue].items()])
                tpc_gene_tissue = nnz_gene_tissue / cell_count_tissue
                gene_tissue_expr_stats = {
                    "n": nnz_gene_tissue,
                    "me": me,
                    "tpc": tpc_gene_tissue,
                }

                result[gene][tissue]["tissue_stats"] = {}
                result[gene][tissue]["tissue_stats"]["aggregated"] = gene_tissue_expr_stats

    return result


def generate_test_inputs_and_expected_outputs(
    genes: list[str],
    organism: str,
    dim_size: int,
    me: float,
    cell_count_per_row_cell_counts_cube: int,
    compare_dim=None,
    cell_ordering_func=forward_cell_type_ordering,
) -> tuple:
    """
    Generates test inputs and expected outputs for the /wmg/v2/query endpoint.

    Arguments
    ---------
    genes: list of gene ontology term IDs

    organism: organism ontology term ID

    dim_size: size of each dimension of the test cube

    me: mean expression value to use for each (gene, tissue, cell_type) combination (scalar)

    cell_count_per_row_cell_counts_cube: num of cells per row in cell_counts cube (scalar)

    compare_dim: dimension to use for compare feature (optional). None if compare isn't used.

    Returns
    -------
    tuple of (request, expected_expression_summary, expected_term_id_labels)

    request: dictionary containing the request body to send to the /wmg/v2/query endpoint

    expected_expression_summary: dictionary containing the expected expression summary values

    expected_term_id_labels: dictionary containing the expected term ID labels
    """
    cell_types = [f"cell_type_ontology_term_id_{i}" for i in range(dim_size)]

    # WMG V2 API does not allow filtering by tissues and therefore the query result
    # includes all tissues
    all_tissues = [f"tissue_ontology_term_id_{i}" for i in range(dim_size)]

    expected_combinations_per_tissue = dim_size ** len(expression_summary_non_indexed_dims)
    cell_count_tissue = cell_count_per_row_cell_counts_cube * expected_combinations_per_tissue

    expected_combinations_per_cell_type = dim_size ** len(
        set(expression_summary_non_indexed_dims).difference({"cell_type_ontology_term_id"})
    )
    nnz_gene_tissue_cell_type = expected_combinations_per_cell_type
    cell_count_tissue_cell_type = expected_combinations_per_cell_type * cell_count_per_row_cell_counts_cube

    compare_terms = []
    cell_counts_tissue_cell_type_compare_dim = 0
    nnz_gene_tissue_cell_type_compare_dim = 0

    if compare_dim:
        compare_terms = [f"{find_dimension_id_from_compare(compare_dim)}_{i}" for i in range(dim_size)]
        cell_counts_tissue_cell_type_compare_dim = cell_count_tissue_cell_type // len(compare_terms)
        nnz_gene_tissue_cell_type_compare_dim = nnz_gene_tissue_cell_type // len(compare_terms)

    expected_term_id_labels = generate_expected_term_id_labels_dictionary(
        genes=genes,
        tissues=all_tissues,
        cell_types=cell_types,
        cell_count_tissue_cell_type=cell_count_tissue_cell_type,
        compare_terms=compare_terms,
        cell_counts_tissue_cell_type_compare_dim=cell_counts_tissue_cell_type_compare_dim,
        cell_ordering_func=cell_ordering_func,
    )
    expected_expression_summary = generate_expected_expression_summary_dictionary(
        genes=genes,
        tissues=all_tissues,
        cell_count_tissue=cell_count_tissue,
        cell_types=cell_types,
        cell_count_tissue_cell_type=cell_count_tissue_cell_type,
        nnz_gene_tissue_cell_type=nnz_gene_tissue_cell_type,
        compare_terms=compare_terms,
        cell_counts_tissue_cell_type_compare_dim=cell_counts_tissue_cell_type_compare_dim,
        nnz_gene_tissue_cell_type_compare_dim=nnz_gene_tissue_cell_type_compare_dim,
        me=me,
    )

    request = dict(filter=dict(gene_ontology_term_ids=genes, organism_ontology_term_id=organism))
    if compare_dim:
        request["compare"] = compare_dim

    return (
        request,
        expected_expression_summary,
        expected_term_id_labels,
    )


# TODO(prathap): Write tests that mock backend.wmg.api.v2.get_dot_plot_data() and
# backend.wmg.api.v2.rollup() so that we can test backend.wmg.api.v2.query() with
# rollup operations.
# see: https://github.com/chanzuckerberg/single-cell-data-portal/issues/4997
class WmgApiV2Tests(unittest.TestCase):
    """
    Tests WMG API endpoints. Tests the flask app only, and not other stack dependencies, such as S3. Builds and uses a
    temporary WMG cube on local filesystem to avoid dependency on localstack S3.
    """

    # TODO(prathap): Write a generalized utility function that extends the
    # comparison capability of sequences to also compare floating point
    # values: https://docs.python.org/3/library/collections.abc.html#module-collections.abc
    def assert_equality_nested_dict_with_floats(self, *, expected, actual, key_path):
        self.assertEqual(expected.keys(), actual.keys(), f"Assertion failure in key path: {key_path}")

        for k, expected_value in expected.items():
            key_path.append(k)

            actual_value = actual[k]

            if isinstance(expected_value, (float, int)) and isinstance(actual_value, (float, int)):
                self.assertEqual(expected_value, approx(actual_value), f"Assertion failure in key path: {key_path}")
            else:
                self.assertIs(type(expected_value), type(actual_value), f"Assertion failure in key path: {key_path}")

                if isinstance(expected_value, dict):
                    self.assert_equality_nested_dict_with_floats(
                        expected=expected_value, actual=actual_value, key_path=key_path
                    )
                else:
                    self.assertEqual(expected_value, actual_value, f"Assertion failure in key path: {key_path}")

            key_path.pop()

        self.assertTrue(True, f"Assertion failure in key path: {key_path}")

    def setUp(self):
        super().setUp()
        with EnvironmentSetup(dict(APP_NAME="corpora-api")):
            self.app = app.test_client(use_cookies=False)

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.maxDiff = None

    @patch("backend.wmg.api.v2.load_snapshot")
    def test__primary_filter_dimensions__returns_200(self, load_snapshot):
        # This test appears to be hitting a TileDB (<=0.13.1) bug and fails (intermittently) if dim_size=3
        with create_temp_wmg_snapshot(dim_size=1) as snapshot:
            # setup up API endpoints to use a mocked cube
            load_snapshot.return_value = snapshot
            response = self.app.get("/wmg/v2/primary_filter_dimensions")

        self.assertEqual(200, response.status_code)

    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
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

            response = self.app.get("/wmg/v2/primary_filter_dimensions")

        expected = dict(
            snapshot_id=test_snapshot_id,
            organism_terms=test_organism_terms,
            tissue_terms=test_tissue_terms,
            gene_terms=test_gene_terms,
        )

        self.assertEqual(expected, json.loads(response.data))

    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
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
            organism = "organism_ontology_term_id_0"

            (request, expected_expression_summary, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, organism, dim_size, 1.0, 10
            )

            response = self.app.post("/wmg/v2/query", json=request)

            self.assertEqual(200, response.status_code)

            expected_response = {
                "snapshot_id": "dummy-snapshot",
                "expression_summary": expected_expression_summary,
                "term_id_labels": expected_term_id_labels,
            }
            self.assert_equality_nested_dict_with_floats(
                expected=expected_response, actual=json.loads(response.data), key_path=[]
            )

    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
    def test__query_no_genes__returns_200_and_correct_response(
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

            # this is the convention used by the FE to indicate no genes are selected
            # see https://github.com/chanzuckerberg/single-cell-data-portal/pull/2729
            genes = ["."]

            organism = "organism_ontology_term_id_0"

            (request, expected_expression_summary, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, organism, dim_size, 1.0, 10
            )

            response = self.app.post("/wmg/v2/query", json=request)

            self.assertEqual(200, response.status_code)

            expected_response = {
                "snapshot_id": "dummy-snapshot",
                "expression_summary": expected_expression_summary,
                "term_id_labels": expected_term_id_labels,
            }
            self.assert_equality_nested_dict_with_floats(
                expected=expected_response, actual=json.loads(response.data), key_path=[]
            )

    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
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
            organism = "organism_ontology_term_id_0"

            (request, expected_expression_summary, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, organism, dim_size, 1.0, 10
            )

            response = self.app.post("/wmg/v2/query", json=request)

            self.assertEqual(200, response.status_code)

            expected = {
                "snapshot_id": "dummy-snapshot",
                "expression_summary": expected_expression_summary,
                "term_id_labels": expected_term_id_labels,
            }

            self.assert_equality_nested_dict_with_floats(
                expected=expected, actual=json.loads(response.data), key_path=[]
            )

    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
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
            organism = "organism_ontology_term_id_0"

            (request, expected_expression_summary, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, organism, dim_size, 1.0, 10, compare_dim="self_reported_ethnicity"
            )

            response = self.app.post("/wmg/v2/query", json=request)

            self.assertEqual(200, response.status_code)

            expected = {
                "snapshot_id": "dummy-snapshot",
                "expression_summary": expected_expression_summary,
                "term_id_labels": expected_term_id_labels,
            }

            self.assert_equality_nested_dict_with_floats(
                expected=expected, actual=json.loads(response.data), key_path=[]
            )

    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
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
            organism = "organism_ontology_term_id_0"

            (request, _, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, organism, dim_size, 1.0, 10, cell_ordering_func=reverse_cell_type_ordering
            )

            response = self.app.post("/wmg/v2/query", json=request)

            self.assertEqual(200, response.status_code)

            expected = expected_term_id_labels["cell_types"]
            self.assertTrue(compare_dicts(expected, json.loads(response.data)["term_id_labels"]["cell_types"]))

    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
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
            organism = "organism_ontology_term_id_0"

            (request, _, expected_term_id_labels) = generate_test_inputs_and_expected_outputs(
                genes, organism, dim_size, 1.0, expected_count, cell_ordering_func=reverse_cell_type_ordering
            )

            response = self.app.post("/wmg/v2/query", json=request)

            self.assertEqual(200, response.status_code)

            expected = expected_term_id_labels["cell_types"]
            self.assertEqual(expected, json.loads(response.data)["term_id_labels"]["cell_types"])

    def test__query_containing_tissue__request_returns_400(self):
        request = dict(
            filter=dict(
                gene_ontology_term_ids=["gene_ontology_term_id_0"],
                organism_ontology_term_id="organism_ontology_term_id_0",
                tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            ),
        )

        response = self.app.post("/wmg/v2/query", json=request)

        self.assertEqual(400, response.status_code)
        self.assertEqual(
            "Additional properties are not allowed ('tissue_ontology_term_ids' was unexpected) - 'filter'",
            json.loads(response.data)["detail"],
        )

    def test__query_empty_request__returns_400(self):
        response = self.app.post("/wmg/v2/query", json={})

        self.assertEqual(400, response.status_code)
        self.assertEqual("'filter' is a required property", json.loads(response.data)["detail"])

    def test__query_missing_gene__request_returns_400(self):
        request = dict(
            filter=dict(
                organism_ontology_term_id="organism_ontology_term_id_0",
            ),
        )

        response = self.app.post("/wmg/v2/query", json=request)

        self.assertEqual(400, response.status_code)
        self.assertEqual(
            "'gene_ontology_term_ids' is a required property - 'filter'", json.loads(response.data)["detail"]
        )

    def test__query_missing_organism__request_returns_400(self):
        request = dict(
            filter=dict(
                gene_ontology_term_ids=["gene_ontology_term_id_0"],
            ),
        )

        response = self.app.post("/wmg/v2/query", json=request)

        self.assertEqual(400, response.status_code)
        self.assertEqual(
            "'organism_ontology_term_id' is a required property - 'filter'", json.loads(response.data)["detail"]
        )

    def test__filter_request_with_query_key_not_specified_in_api_spec_returns_400(self):
        filter_dict = dict(
            # these don't matter for the expected result
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            # these matter for the expected result
            development_stage_ontology_term_ids=["development_stage_ontology_term_id_0"],
            self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_0"],
        )

        filter_request = dict(
            filter=filter_dict,
            is_rollup=True,
        )

        response = self.app.post("/wmg/v2/filters", json=filter_request)

        self.assertEqual(400, response.status_code)
        self.assertEqual(
            "Additional properties are not allowed ('is_rollup' was unexpected)", json.loads(response.data)["detail"]
        )

    @patch("backend.wmg.api.v2.fetch_datasets_metadata")
    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
    def test__filter_request_with_filter_dims__returns_valid_filter_dims__base_case(
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

            filter_0_request = dict(filter=filter_0)

            response = self.app.post("/wmg/v2/filters", json=filter_0_request)

            expected_filters = {
                "cell_type_terms": [{"cell_type_ontology_term_id_0": "cell_type_ontology_term_id_0_label"}],
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
                "publication_citations": [],
                "self_reported_ethnicity_terms": [
                    {"self_reported_ethnicity_ontology_term_id_0": "self_reported_ethnicity_ontology_term_id_0_label"}
                ],
                "sex_terms": [],
                "tissue_terms": [{"tissue_ontology_term_id_0": "tissue_ontology_term_id_0_label"}],
            }
            self.assertEqual(json.loads(response.data)["filter_dims"], expected_filters)

    @patch("backend.wmg.api.v2.fetch_datasets_metadata")
    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
    def test__filter_request_with_filter_dims__returns_valid_filter_dims(
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

                filter_0_request = dict(filter=filter_0)

                filter_0_no_dev_stage_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=[],
                    self_reported_ethnicity_ontology_term_ids=["self_reported_ethnicity_ontology_term_id_0"],
                )
                filter_0_no_dev_stage_request = dict(filter=filter_0_no_dev_stage_filter)

                filter_0_no_ethnicity_filter = dict(
                    # these don't matter for the expected result
                    gene_ontology_term_ids=["gene_ontology_term_id_0"],
                    organism_ontology_term_id="organism_ontology_term_id_0",
                    tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
                    # these matter for the expected result
                    development_stage_ontology_term_ids=["development_stage_ontology_term_id_0"],
                    self_reported_ethnicity_ontology_term_ids=[],
                )
                filter_0_no_ethnicity_request = dict(filter=filter_0_no_ethnicity_filter)
                # the values for dev_stage terms when a dev stage filter is included should match the values returned
                # if no filter is passed in for dev stage
                response = self.app.post("/wmg/v2/filters", json=filter_0_request)
                dev_stage_terms = json.loads(response.data)["filter_dims"]["development_stage_terms"]
                self_reported_ethnicity_terms = json.loads(response.data)["filter_dims"][
                    "self_reported_ethnicity_terms"
                ]

                no_dev_stage_filter_response = self.app.post("/wmg/v2/filters", json=filter_0_no_dev_stage_request)
                dev_stage_terms_if_no_dev_stage_filters = json.loads(no_dev_stage_filter_response.data)["filter_dims"][
                    "development_stage_terms"
                ]

                no_ethnicity_filter_response = self.app.post("/wmg/v2/filters", json=filter_0_no_ethnicity_request)
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
                self_reported_ethnicity_0_request = dict(filter=self_reported_ethnicity_0_filter)
                response = self.app.post("/wmg/v2/filters", json=self_reported_ethnicity_0_request)
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
                self_reported_ethnicity_1_request = dict(filter=self_reported_ethnicity_1_filter)
                response = self.app.post("/wmg/v2/filters", json=self_reported_ethnicity_1_request)
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
                self_reported_ethnicity_1_dev_1_request = dict(filter=self_reported_ethnicity_1_dev_1_filter)

                response = self.app.post("/wmg/v2/filters", json=self_reported_ethnicity_1_dev_1_request)
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
                self_reported_ethnicity_2_request = dict(filter=self_reported_ethnicity_2_filter)
                eth_2_dev_2_request = dict(filter=self_reported_ethnicity_2_dev_2_filter)
                response = self.app.post("/wmg/v2/filters", json=self_reported_ethnicity_2_request)
                dev_stage_terms_eth_2_no_dev_filter = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["development_stage_terms"]
                )
                self_reported_ethnicity_terms_eth_2_no_dev_filter = _sort_list_by_dictionary_keys(
                    json.loads(response.data)["filter_dims"]["self_reported_ethnicity_terms"]
                )
                self.assertEqual(expected_development_stage_terms, dev_stage_terms_eth_2_no_dev_filter)
                self.assertEqual(all_self_reported_ethnicity_terms, self_reported_ethnicity_terms_eth_2_no_dev_filter)

                response = self.app.post("/wmg/v2/filters", json=eth_2_dev_2_request)
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

    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
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

            response = self.app.post("/wmg/v2/markers", json=request)
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

    @patch("backend.wmg.api.v2.gene_term_label")
    @patch("backend.wmg.api.v2.ontology_term_label")
    @patch("backend.wmg.api.v2.load_snapshot")
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

            response = self.app.post("/wmg/v2/markers", json=request)
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
