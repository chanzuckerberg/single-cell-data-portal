import unittest

from backend.wmg.data.query import WmgQueryCriteria, WmgQuery
from backend.wmg.data.schema import cube_non_indexed_dims
from unit.backend.wmg.fixtures.cube import create_temp_cube, all_ones_attr_values


class QueryTest(unittest.TestCase):
    def test__query_single_gene_organism_and_tissue__returns_correct_result(self):
        criteria = WmgQueryCriteria(
            organism_term_id="organism_ontology_term_id_0", tissue_term_ids=["tissue_ontology_term_id_0"]
        )

        dim_size = 3
        with create_temp_cube(dim_size=dim_size, attr_vals_fn=all_ones_attr_values) as cube:
            result = WmgQuery(cube).execute(criteria)

        # sanity check the expected values of the `n` stat for each data viz point; if this fails, the cube test
        # fixture may have changed (e.g. TileDB Array schema) or the logic for creating the test cube fixture has
        # changed
        expected_cell_count_per_cell_type = dim_size ** len(
            set(cube_non_indexed_dims).difference({"cell_type_ontology_term_id"})
        )
        assert expected_cell_count_per_cell_type == 729

        expected = {
            "n_cells": {
                (b"gene_term_id_0", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_0"): 729,
                (b"gene_term_id_0", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_1"): 729,
                (b"gene_term_id_0", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_2"): 729,
                (b"gene_term_id_1", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_0"): 729,
                (b"gene_term_id_1", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_1"): 729,
                (b"gene_term_id_1", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_2"): 729,
                (b"gene_term_id_2", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_0"): 729,
                (b"gene_term_id_2", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_1"): 729,
                (b"gene_term_id_2", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_2"): 729,
            },
            "nnz": {
                (b"gene_term_id_0", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_0"): 729,
                (b"gene_term_id_0", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_1"): 729,
                (b"gene_term_id_0", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_2"): 729,
                (b"gene_term_id_1", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_0"): 729,
                (b"gene_term_id_1", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_1"): 729,
                (b"gene_term_id_1", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_2"): 729,
                (b"gene_term_id_2", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_0"): 729,
                (b"gene_term_id_2", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_1"): 729,
                (b"gene_term_id_2", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_2"): 729,
            },
            "sum": {
                (b"gene_term_id_0", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_0"): 729.0,
                (b"gene_term_id_0", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_1"): 729.0,
                (b"gene_term_id_0", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_2"): 729.0,
                (b"gene_term_id_1", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_0"): 729.0,
                (b"gene_term_id_1", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_1"): 729.0,
                (b"gene_term_id_1", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_2"): 729.0,
                (b"gene_term_id_2", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_0"): 729.0,
                (b"gene_term_id_2", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_1"): 729.0,
                (b"gene_term_id_2", b"tissue_ontology_term_id_0", b"cell_type_ontology_term_id_2"): 729.0,
            },
        }

        self.assertEqual(expected, result.to_dict())
