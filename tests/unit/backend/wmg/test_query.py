import unittest
from typing import NamedTuple

from backend.wmg.api.v1 import agg_cell_type_counts, agg_tissue_counts, get_dot_plot_data
from backend.wmg.data.query import (
    FmgQueryCriteria,
    MarkerGeneQueryCriteria,
    WmgQuery,
    WmgQueryCriteria,
    retrieve_top_n_markers,
)
from backend.wmg.data.schemas.cube_schema import expression_summary_non_indexed_dims
from tests.unit.backend.wmg.fixtures.test_snapshot import (
    all_ones_expression_summary_values,
    all_tens_cell_counts_values,
    all_X_cell_counts_values,
    create_temp_wmg_snapshot,
    load_realistic_test_snapshot,
)

TEST_SNAPSHOT = "realistic-test-snapshot"

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
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            result = q.marker_genes(criteria)
            marker_genes = retrieve_top_n_markers(result, "ttest", 10)
            expected = [
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
            ]
            self.assertEqual(marker_genes, expected)

    def test__query_marker_genes_cube__returns_correct_all_markers(self):
        criteria = MarkerGeneQueryCriteria(
            tissue_ontology_term_id="UBERON:0002048",
            cell_type_ontology_term_id="CL:0000786",
            organism_ontology_term_id="NCBITaxon:9606",
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            result = q.marker_genes(criteria)
            marker_genes = retrieve_top_n_markers(result, "ttest", 0)
            expected = [
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
                {"gene_ontology_term_id": "ENSG00000173110", "p_value": 0.0, "effect_size": 0.388088583946228},
                {"gene_ontology_term_id": "ENSG00000171067", "p_value": 0.0, "effect_size": 0.19824208319187164},
                {"gene_ontology_term_id": "ENSG00000188372", "p_value": 0.0, "effect_size": 0.12499640136957169},
                {"gene_ontology_term_id": "ENSG00000103227", "p_value": 0.0, "effect_size": 0.11911825835704803},
                {"gene_ontology_term_id": "ENSG00000133460", "p_value": 0.0, "effect_size": 0.08775933086872101},
                {"gene_ontology_term_id": "ENSG00000171311", "p_value": 0.0, "effect_size": 0.06580302119255066},
                {"gene_ontology_term_id": "ENSG00000051108", "p_value": 0.0, "effect_size": 0.06569653749465942},
                {"gene_ontology_term_id": "ENSG00000170542", "p_value": 0.0, "effect_size": 0.06567766517400742},
                {"gene_ontology_term_id": "ENSG00000031698", "p_value": 0.0, "effect_size": 0.033098865300416946},
                {"gene_ontology_term_id": "ENSG00000102158", "p_value": 0.0, "effect_size": -0.004498007707297802},
                {"gene_ontology_term_id": "ENSG00000245060", "p_value": 0.0, "effect_size": -0.004812875296920538},
                {"gene_ontology_term_id": "ENSG00000148690", "p_value": 0.0, "effect_size": -0.004902553744614124},
                {"gene_ontology_term_id": "ENSG00000163568", "p_value": 0.0, "effect_size": -0.0326923169195652},
                {"gene_ontology_term_id": "ENSG00000244462", "p_value": 0.0, "effect_size": -0.06741894781589508},
                {"gene_ontology_term_id": "ENSG00000182162", "p_value": 0.0, "effect_size": -0.08986986428499222},
                {"gene_ontology_term_id": "ENSG00000116127", "p_value": 0.0, "effect_size": -0.09243863821029663},
                {"gene_ontology_term_id": "ENSG00000102007", "p_value": 0.0, "effect_size": -0.110072061419487},
                {"gene_ontology_term_id": "ENSG00000117133", "p_value": 0.0, "effect_size": -0.13198323547840118},
                {"gene_ontology_term_id": "ENSG00000084072", "p_value": 0.0, "effect_size": -0.13892346620559692},
                {"gene_ontology_term_id": "ENSG00000109790", "p_value": 0.0, "effect_size": -0.14885607361793518},
                {"gene_ontology_term_id": "ENSG00000127526", "p_value": 0.0, "effect_size": -0.1617552489042282},
                {"gene_ontology_term_id": "ENSG00000092820", "p_value": 0.0, "effect_size": -0.1630970537662506},
                {"gene_ontology_term_id": "ENSG00000184661", "p_value": 0.0, "effect_size": -0.1941165179014206},
                {"gene_ontology_term_id": "ENSG00000164896", "p_value": 0.0, "effect_size": -0.2109561562538147},
                {"gene_ontology_term_id": "ENSG00000125740", "p_value": 0.0, "effect_size": -0.2233702838420868},
                {"gene_ontology_term_id": "ENSG00000140553", "p_value": 0.0, "effect_size": -0.2270195335149765},
                {"gene_ontology_term_id": "ENSG00000145781", "p_value": 0.0, "effect_size": -0.2825307548046112},
                {"gene_ontology_term_id": "ENSG00000134970", "p_value": 0.0, "effect_size": -0.33043113350868225},
                {"gene_ontology_term_id": "ENSG00000109171", "p_value": 0.0, "effect_size": -0.35054484009742737},
                {"gene_ontology_term_id": "ENSG00000132254", "p_value": 0.0, "effect_size": -0.3540075719356537},
                {"gene_ontology_term_id": "ENSG00000143442", "p_value": 0.0, "effect_size": -0.35437655448913574},
                {"gene_ontology_term_id": "ENSG00000148824", "p_value": 0.0, "effect_size": -0.36991795897483826},
                {"gene_ontology_term_id": "ENSG00000163029", "p_value": 0.0, "effect_size": -0.38332968950271606},
                {"gene_ontology_term_id": "ENSG00000116473", "p_value": 0.0, "effect_size": -0.40905269980430603},
                {"gene_ontology_term_id": "ENSG00000002330", "p_value": 0.0, "effect_size": -0.4481770396232605},
                {"gene_ontology_term_id": "ENSG00000186818", "p_value": 0.0, "effect_size": -0.4491439163684845},
                {"gene_ontology_term_id": "ENSG00000152061", "p_value": 0.0, "effect_size": -0.5008719563484192},
                {"gene_ontology_term_id": "ENSG00000173992", "p_value": 0.0, "effect_size": -0.5449978113174438},
                {"gene_ontology_term_id": "ENSG00000152818", "p_value": 0.0, "effect_size": -0.5866119861602783},
                {"gene_ontology_term_id": "ENSG00000265206", "p_value": 0.0, "effect_size": -0.5957227349281311},
                {"gene_ontology_term_id": "ENSG00000122483", "p_value": 0.0, "effect_size": -0.6621347069740295},
                {"gene_ontology_term_id": "ENSG00000070831", "p_value": 0.0, "effect_size": -0.6757633090019226},
                {"gene_ontology_term_id": "ENSG00000159128", "p_value": 0.0, "effect_size": -0.7019357085227966},
                {"gene_ontology_term_id": "ENSG00000186184", "p_value": 0.0, "effect_size": -0.7033576369285583},
                {"gene_ontology_term_id": "ENSG00000106367", "p_value": 0.0, "effect_size": -0.7177350521087646},
                {"gene_ontology_term_id": "ENSG00000106554", "p_value": 0.0, "effect_size": -0.7232860326766968},
                {"gene_ontology_term_id": "ENSG00000100445", "p_value": 0.0, "effect_size": -0.7570595145225525},
                {"gene_ontology_term_id": "ENSG00000135940", "p_value": 0.0, "effect_size": -0.7782347202301025},
                {"gene_ontology_term_id": "ENSG00000119326", "p_value": 0.0, "effect_size": -0.8128666877746582},
                {"gene_ontology_term_id": "ENSG00000122417", "p_value": 0.0, "effect_size": -0.8753214478492737},
                {"gene_ontology_term_id": "ENSG00000168028", "p_value": 0.0, "effect_size": -0.9302230477333069},
                {"gene_ontology_term_id": "ENSG00000161547", "p_value": 0.0, "effect_size": -1.0066393613815308},
                {"gene_ontology_term_id": "ENSG00000244187", "p_value": 0.0, "effect_size": -1.0319607257843018},
                {"gene_ontology_term_id": "ENSG00000113580", "p_value": 0.0, "effect_size": -1.0370265245437622},
                {"gene_ontology_term_id": "ENSG00000075415", "p_value": 0.0, "effect_size": -1.2570115327835083},
                {"gene_ontology_term_id": "ENSG00000276043", "p_value": 0.0, "effect_size": -1.298696756362915},
                {"gene_ontology_term_id": "ENSG00000133112", "p_value": 0.0, "effect_size": -1.3088808059692383},
                {"gene_ontology_term_id": "ENSG00000170296", "p_value": 0.0, "effect_size": -1.9430065155029297},
            ]
            self.assertEqual(marker_genes, expected)

    def test__query_expression_summary_fmg_cube__returns_correct_results(self):
        criteria = FmgQueryCriteria(
            gene_ontology_term_ids=["ENSG00000238042", "ENSG00000168028"],
            organism_ontology_term_id="NCBITaxon:9606",
            tissue_ontology_term_ids=["UBERON:0002048"],
            cell_type_ontology_term_ids=["CL:0000786"],
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            query_result = q.expression_summary_fmg(criteria)
            query_sum = list(query_result[["sum", "sqsum", "nnz", "nnz_thr"]].sum())
            expected = [85958.453125, 241523.03125, 37404.0, 36434.0]
            [self.assertAlmostEqual(query_sum[i], expected[i], places=3) for i in range(len(query_sum))]

    def test__query_expression_summary_default_cube__returns_correct_results(self):
        criteria = WmgQueryCriteria(
            gene_ontology_term_ids=["ENSG00000238042", "ENSG00000168028"],
            organism_ontology_term_id="NCBITaxon:9606",
            tissue_ontology_term_ids=["UBERON:0002048"],
        )
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            q = WmgQuery(snapshot)
            query_result = q.expression_summary_default(criteria)
            query_sum = list(query_result[["sum", "nnz"]].sum())
            expected = [553378.625, 261191.0]
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
