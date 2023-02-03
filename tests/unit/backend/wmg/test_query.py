import unittest
from typing import NamedTuple

from backend.wmg.api.v1 import get_dot_plot_data, agg_cell_type_counts, agg_tissue_counts
from backend.wmg.data.query import (
    WmgQueryCriteria,
    WmgQuery,
    MarkerGeneQueryCriteria,
    FmgQueryCriteria,
    retrieve_top_n_markers,
)
from backend.wmg.data.schemas.cube_schema import expression_summary_non_indexed_dims
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    create_temp_wmg_snapshot,
    load_test_fmg_snapshot,
    all_ones_expression_summary_values,
    all_tens_cell_counts_values,
    all_X_cell_counts_values,
)

TEST_SNAPSHOT = "test-fmg-snapshot"

ALL_INDEXED_DIMS_FOR_QUERY = [
    "gene_ontology_term_ids",
    "tissue_ontology_term_ids",
    "tissue_original_ontology_term_ids",
    "organism_ontology_term_id",
]

# TODO: Test build_* methods separately in test_v1.py.  This package's unit tests need only test the raw results of
#  query methods


class QueryTest(unittest.TestCase):
    def test__query_marker_genes_cube__returns_correct_top_10_markers(self):
        criteria = MarkerGeneQueryCriteria(
            tissue_ontology_term_id="UBERON:0002048",
            cell_type_ontology_term_id="CL:0000786",
            organism_ontology_term_id="NCBITaxon:9606",
        )
        with load_test_fmg_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            result = q.marker_genes(criteria)
            marker_genes = retrieve_top_n_markers(result, "ttest", 10)
            expected = [
                {"gene_ontology_term_id": "ENSG00000132465", "p_value": 0.0, "effect_size": 2.661743640899658},
                {"gene_ontology_term_id": "ENSG00000180879", "p_value": 0.0, "effect_size": 2.1092216968536377},
                {"gene_ontology_term_id": "ENSG00000170476", "p_value": 0.0, "effect_size": 2.072816848754883},
                {"gene_ontology_term_id": "ENSG00000134285", "p_value": 0.0, "effect_size": 1.8734184503555298},
                {"gene_ontology_term_id": "ENSG00000166562", "p_value": 0.0, "effect_size": 1.623068928718567},
                {"gene_ontology_term_id": "ENSG00000099958", "p_value": 0.0, "effect_size": 1.5588932037353516},
                {"gene_ontology_term_id": "ENSG00000211592", "p_value": 0.0, "effect_size": 1.4483562707901},
                {"gene_ontology_term_id": "ENSG00000051108", "p_value": 0.0, "effect_size": 1.2816050052642822},
                {"gene_ontology_term_id": "ENSG00000118363", "p_value": 0.0, "effect_size": 1.0649874210357666},
                {"gene_ontology_term_id": "ENSG00000100219", "p_value": 0.0, "effect_size": 0.46257320046424866},
            ]

            self.assertEqual(marker_genes, expected)

    def test__query_marker_genes_cube__returns_correct_all_markers(self):
        criteria = MarkerGeneQueryCriteria(
            tissue_ontology_term_id="UBERON:0002048",
            cell_type_ontology_term_id="CL:0000786",
            organism_ontology_term_id="NCBITaxon:9606",
        )
        with load_test_fmg_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            result = q.marker_genes(criteria)
            marker_genes = retrieve_top_n_markers(result, "ttest", 0)
            expected = [
                {"gene_ontology_term_id": "ENSG00000132465", "p_value": 0.0, "effect_size": 2.661743640899658},
                {"gene_ontology_term_id": "ENSG00000180879", "p_value": 0.0, "effect_size": 2.1092216968536377},
                {"gene_ontology_term_id": "ENSG00000170476", "p_value": 0.0, "effect_size": 2.072816848754883},
                {"gene_ontology_term_id": "ENSG00000134285", "p_value": 0.0, "effect_size": 1.8734184503555298},
                {"gene_ontology_term_id": "ENSG00000166562", "p_value": 0.0, "effect_size": 1.623068928718567},
                {"gene_ontology_term_id": "ENSG00000099958", "p_value": 0.0, "effect_size": 1.5588932037353516},
                {"gene_ontology_term_id": "ENSG00000211592", "p_value": 0.0, "effect_size": 1.4483562707901},
                {"gene_ontology_term_id": "ENSG00000051108", "p_value": 0.0, "effect_size": 1.2816050052642822},
                {"gene_ontology_term_id": "ENSG00000118363", "p_value": 0.0, "effect_size": 1.0649874210357666},
                {"gene_ontology_term_id": "ENSG00000100219", "p_value": 0.0, "effect_size": 0.46257320046424866},
                {"gene_ontology_term_id": "ENSG00000173110", "p_value": 0.0, "effect_size": 0.04467674344778061},
                {"gene_ontology_term_id": "ENSG00000125844", "p_value": 0.0, "effect_size": -0.018479663878679276},
                {"gene_ontology_term_id": "ENSG00000102007", "p_value": 0.0, "effect_size": -0.03636987507343292},
                {"gene_ontology_term_id": "ENSG00000031698", "p_value": 0.0, "effect_size": -0.10610716789960861},
                {"gene_ontology_term_id": "ENSG00000102158", "p_value": 0.0, "effect_size": -0.13100826740264893},
                {"gene_ontology_term_id": "ENSG00000186818", "p_value": 0.0, "effect_size": -0.18612566590309143},
                {"gene_ontology_term_id": "ENSG00000159128", "p_value": 0.0, "effect_size": -0.5121894478797913},
                {"gene_ontology_term_id": "ENSG00000125740", "p_value": 0.0, "effect_size": -0.6192882657051086},
                {"gene_ontology_term_id": "ENSG00000116473", "p_value": 0.0, "effect_size": -0.6428257822990417},
                {"gene_ontology_term_id": "ENSG00000133112", "p_value": 0.0, "effect_size": -0.6987436413764954},
                {"gene_ontology_term_id": "ENSG00000092820", "p_value": 0.0, "effect_size": -0.7473370432853699},
                {"gene_ontology_term_id": "ENSG00000161547", "p_value": 0.0, "effect_size": -0.776445746421814},
                {"gene_ontology_term_id": "ENSG00000186184", "p_value": 0.0, "effect_size": -0.7776976823806763},
                {"gene_ontology_term_id": "ENSG00000070831", "p_value": 0.0, "effect_size": -0.8825482130050659},
                {"gene_ontology_term_id": "ENSG00000135940", "p_value": 0.0, "effect_size": -0.8879522085189819},
                {"gene_ontology_term_id": "ENSG00000168028", "p_value": 0.0, "effect_size": -1.3053345680236816},
                {"gene_ontology_term_id": "ENSG00000075415", "p_value": 0.0, "effect_size": -1.390924096107483},
                {"gene_ontology_term_id": "ENSG00000170296", "p_value": 0.0, "effect_size": -1.5161384344100952},
            ]
            self.assertEqual(marker_genes, expected)

    def test__query_expression_summary_fmg_cube__returns_correct_results(self):
        criteria = FmgQueryCriteria(
            gene_ontology_term_ids=["ENSG00000238042", "ENSG00000168028"],
            organism_ontology_term_id="NCBITaxon:9606",
            tissue_ontology_term_ids=["UBERON:0002048"],
            cell_type_ontology_term_ids=["CL:0000786"],
        )
        with load_test_fmg_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            query_result = q.expression_summary_fmg(criteria)
            query_sum = list(query_result[["sum", "sqsum", "nnz", "nnz_thr"]].sum())
            expected = [28538.255859375, 85875.046875, 11312.0, 11185.0]
            [self.assertAlmostEqual(query_sum[i], expected[i], places=3) for i in range(len(query_sum))]

    def test__query_all_indexed_dims_single_value__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_1",
            tissue_ontology_term_ids=["tissue_ontology_term_id_2"],
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            q = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(q.expression_summary(criteria), q.cell_counts(criteria))

        # sanity check the expected value of the stats (nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        not_used_cube_indexed_dims = [0 if criteria.dict()[dim_name] else 1 for dim_name in ALL_INDEXED_DIMS_FOR_QUERY]

        expected_cell_count_per_cell_type = dim_size ** (
            len(expression_summary_non_indexed_dims) + sum(not_used_cube_indexed_dims) - 1
        )

        expected_cell_count_per_tissue = 10 * (
            dim_size ** (len(expression_summary_non_indexed_dims) + sum(not_used_cube_indexed_dims))
        )

        assert expected_cell_count_per_cell_type == 2187
        assert expected_cell_count_per_tissue == 65610

        expected = [
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
        ]

        self.assertEqual(
            expected,
            sorted(
                result.to_dict("records"),
                key=lambda r: (
                    r["gene_ontology_term_id"],
                    r["tissue_ontology_term_id"],
                    r["cell_type_ontology_term_id"],
                ),
            ),
        )

    def test__query_all_indexed_dims_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0", "gene_ontology_term_id_2"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_1", "tissue_ontology_term_id_2"],
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            q = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(q.expression_summary(criteria), q.cell_counts(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        not_used_cube_indexed_dims = [0 if criteria.dict()[dim_name] else 1 for dim_name in ALL_INDEXED_DIMS_FOR_QUERY]
        expected_cell_count_per_cell_type = dim_size ** (len(expression_summary_non_indexed_dims) - 1 + 1)
        expected_cell_count_per_tissue = 10 * (
            dim_size ** (len(expression_summary_non_indexed_dims) + sum(not_used_cube_indexed_dims))
        )

        assert expected_cell_count_per_cell_type == 2187
        assert expected_cell_count_per_tissue == 65610

        expected = [
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_1",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_2",
                "tissue_ontology_term_id": "tissue_ontology_term_id_2",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells_cell_type": 21870,
                "n_cells_tissue": 65610,
                "nnz": 2187,
                "sum": 2187.0,
            },
        ]

        self.assertEqual(
            expected,
            sorted(
                result.to_dict("records"),
                key=lambda r: (
                    r["gene_ontology_term_id"],
                    r["tissue_ontology_term_id"],
                    r["cell_type_ontology_term_id"],
                ),
            ),
        )

    def test__query_agg_cell_type_counts__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_1",
            tissue_ontology_term_ids=[
                "tissue_ontology_term_id_0",
                "tissue_ontology_term_id_1",
                "tissue_ontology_term_id_2",
            ],
        )

        dim_size = 3
        expected_count = 42
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=lambda coords: all_X_cell_counts_values(coords, expected_count),
        ) as snapshot:
            q = WmgQuery(snapshot)
            result = agg_cell_type_counts(q.cell_counts(criteria))

        not_used_cube_indexed_dims = [0 if criteria.dict()[dim_name] else 1 for dim_name in ALL_INDEXED_DIMS_FOR_QUERY]
        expected_n_combinations = dim_size ** (
            len(expression_summary_non_indexed_dims) + sum(not_used_cube_indexed_dims) - 1
        )

        assert expected_n_combinations == 2187

        # after aggregating, we will get three tissues, and three cell types per tissue,
        # with 729 * expected_count total cells
        expected = [{"n_cells_cell_type": 2187 * expected_count}] * len(criteria.tissue_ontology_term_ids) * dim_size
        self.assertEqual(expected, result.to_dict("records"))

    def test__query_agg_tissue_counts__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_1",
            tissue_ontology_term_ids=[
                "tissue_ontology_term_id_0",
                "tissue_ontology_term_id_1",
                "tissue_ontology_term_id_2",
            ],
        )

        dim_size = 3
        expected_count = 42
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=lambda coords: all_X_cell_counts_values(coords, expected_count),
        ) as snapshot:
            q = WmgQuery(snapshot)
            result = agg_tissue_counts(q.cell_counts(criteria))

        not_used_cube_indexed_dims = [0 if criteria.dict()[dim_name] else 1 for dim_name in ALL_INDEXED_DIMS_FOR_QUERY]
        expected_n_combinations = dim_size ** (
            len(expression_summary_non_indexed_dims) + sum(not_used_cube_indexed_dims) - 1
        )

        assert expected_n_combinations == 2187

        # after aggregating, we will get three tissues,
        # with 729 * expected_count * (# cell types per tissue = 3) total cells
        expected = [
            {
                "n_cells_tissue": 2187 * expected_count * dim_size,
            }
        ] * len(criteria.tissue_ontology_term_ids)
        self.assertEqual(expected, result.to_dict("records"))

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.maxDiff = None

    def test__query_non_indexed_dim_single_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            dataset_ids=["dataset_id_1"],  # <-- non-indexed dim, single-valued
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            q = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(q.expression_summary(criteria), q.cell_counts(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed

        not_used_cube_indexed_dims = [0 if criteria.dict()[dim_name] else 1 for dim_name in ALL_INDEXED_DIMS_FOR_QUERY]
        expected_cell_count_per_cell_type = dim_size ** (
            len(expression_summary_non_indexed_dims) - 2 + sum(not_used_cube_indexed_dims)
        )
        expected_cell_count_per_tissue = 10 * (
            dim_size ** (len(expression_summary_non_indexed_dims) - 1 + sum(not_used_cube_indexed_dims))
        )

        assert expected_cell_count_per_cell_type == 729
        assert expected_cell_count_per_tissue == 21870

        expected = [
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells_cell_type": 7290,
                "n_cells_tissue": 21870,
                "nnz": 729,
                "sum": 729.0,
            },
        ]

        self.assertEqual(
            expected,
            sorted(
                result.to_dict("records"),
                key=lambda r: (
                    r["gene_ontology_term_id"],
                    r["tissue_ontology_term_id"],
                    r["cell_type_ontology_term_id"],
                ),
            ),
        )

    def test__query_non_indexed_dim_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            dataset_ids=["dataset_id_1", "dataset_id_0"],  # <-- non-indexed dim, multi-valued
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            q = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(q.expression_summary(criteria), q.cell_counts(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        not_used_cube_indexed_dims = [0 if criteria.dict()[dim_name] else 1 for dim_name in ALL_INDEXED_DIMS_FOR_QUERY]
        expected_cell_count_per_cell_type = (
            dim_size ** (len(expression_summary_non_indexed_dims) - 2 + sum(not_used_cube_indexed_dims)) * 2
        )

        expected_cell_count_per_tissue = 10 * (
            dim_size ** (len(expression_summary_non_indexed_dims) - 1 + sum(not_used_cube_indexed_dims)) * 2
        )

        assert expected_cell_count_per_cell_type == 1458
        assert expected_cell_count_per_tissue == 43740

        expected = [
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells_cell_type": 14580,
                "n_cells_tissue": 43740,
                "nnz": 1458,
                "sum": 1458.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells_cell_type": 14580,
                "n_cells_tissue": 43740,
                "nnz": 1458,
                "sum": 1458.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells_cell_type": 14580,
                "n_cells_tissue": 43740,
                "nnz": 1458,
                "sum": 1458.0,
            },
        ]

        self.assertEqual(
            expected,
            sorted(
                result.to_dict("records"),
                key=lambda r: (
                    r["gene_ontology_term_id"],
                    r["tissue_ontology_term_id"],
                    r["cell_type_ontology_term_id"],
                ),
            ),
        )

    def test__query_non_indexed_dim_single_and_multi_valued__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["gene_ontology_term_id_0"],
            organism_ontology_term_id="organism_ontology_term_id_0",
            tissue_ontology_term_ids=["tissue_ontology_term_id_0"],
            self_reported_ethnicity_ontology_term_ids=[
                "self_reported_ethnicity_ontology_term_id_1"
            ],  # <-- non-indexed dim, single-valued
            dataset_ids=["dataset_id_1", "dataset_id_0"],  # <-- non-indexed dim, multi-valued
        )

        dim_size = 3
        with create_temp_wmg_snapshot(
            dim_size=dim_size,
            expression_summary_vals_fn=all_ones_expression_summary_values,
            cell_counts_generator_fn=all_tens_cell_counts_values,
        ) as snapshot:
            q = WmgQuery(snapshot)
            result, _ = get_dot_plot_data(q.expression_summary(criteria), q.cell_counts(criteria))

        # sanity check the expected value of the stats (n_cells, nnz, sum) for each data viz point; if this fails, the
        # cube test fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture
        # has changed
        not_used_cube_indexed_dims = [0 if criteria.dict()[dim_name] else 1 for dim_name in ALL_INDEXED_DIMS_FOR_QUERY]
        expected_cell_count_per_cell_type = (
            dim_size ** (len(expression_summary_non_indexed_dims) - 3 + sum(not_used_cube_indexed_dims)) * 1 * 2
        )
        expected_cell_count_per_tissue = 10 * (
            dim_size ** (len(expression_summary_non_indexed_dims) - 2 + sum(not_used_cube_indexed_dims)) * 1 * 2
        )

        assert expected_cell_count_per_cell_type == 486
        assert expected_cell_count_per_tissue == 14580

        expected = [
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_0",
                "n_cells_cell_type": 4860,
                "n_cells_tissue": 14580,
                "nnz": 486,
                "sum": 486.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_1",
                "n_cells_cell_type": 4860,
                "n_cells_tissue": 14580,
                "nnz": 486,
                "sum": 486.0,
            },
            {
                "gene_ontology_term_id": "gene_ontology_term_id_0",
                "tissue_ontology_term_id": "tissue_ontology_term_id_0",
                "cell_type_ontology_term_id": "cell_type_ontology_term_id_2",
                "n_cells_cell_type": 4860,
                "n_cells_tissue": 14580,
                "nnz": 486,
                "sum": 486.0,
            },
        ]

        self.assertEqual(
            expected,
            sorted(
                result.to_dict("records"),
                key=lambda r: (
                    r["gene_ontology_term_id"],
                    r["tissue_ontology_term_id"],
                    r["cell_type_ontology_term_id"],
                ),
            ),
        )


class QueryPrimaryFilterDimensionsTest(unittest.TestCase):
    def test__single_dimension__returns_all_dimension_and_terms(self):
        dim_size = 3
        with create_temp_wmg_snapshot(dim_size=dim_size) as snapshot:
            q = WmgQuery(snapshot)
            result = q.list_primary_filter_dimension_term_ids("tissue_ontology_term_id")
            self.assertEquals(
                ["tissue_ontology_term_id_0", "tissue_ontology_term_id_1", "tissue_ontology_term_id_2"], result
            )

    def test__multiple_dimensions__returns_all_dimensions_and_terms_as_tuples(self):
        dim_size = 3

        def exclude_one_tissue_per_organism(logical_coord: NamedTuple) -> bool:
            return logical_coord.tissue_ontology_term_id == logical_coord.organism_ontology_term_id.replace(
                "organism", "tissue"
            )

        with create_temp_wmg_snapshot(
            dim_size=dim_size, exclude_logical_coord_fn=exclude_one_tissue_per_organism
        ) as snapshot:
            q = WmgQuery(snapshot)
            result = q.list_grouped_primary_filter_dimensions_term_ids(
                "tissue_ontology_term_id", "organism_ontology_term_id"
            )
            self.assertEquals(
                {
                    "organism_ontology_term_id_0": ["tissue_ontology_term_id_1", "tissue_ontology_term_id_2"],
                    "organism_ontology_term_id_1": ["tissue_ontology_term_id_0", "tissue_ontology_term_id_2"],
                    "organism_ontology_term_id_2": ["tissue_ontology_term_id_0", "tissue_ontology_term_id_1"],
                },
                result,
            )
