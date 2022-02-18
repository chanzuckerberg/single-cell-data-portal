import unittest

import numpy as np

from unit.backend.wmg.fixtures.cube import create_temp_cube


class QueryTest(unittest.TestCase):

    def test__ad_hoc_query(self):
        with create_temp_cube() as cube:
            feature_ids = cube.query(dims=["feature_id"],
                                     attrs=['n_cells'],
                                     use_arrow=False).df[:]
            print(feature_ids)

    def test__query_single_gene_organism_and_tissue__returns_correct_result(self):
        def attr_values(coords):
            attr_vals = {
                'nnz': np.zeros(len(coords)),
                'n_cells': np.zeros(len(coords)),
                'sum': np.zeros(len(coords))
            }
            for i, coord in enumerate(coords):
                if (coord[0] == 'feature_id_0'
                    and coord[2] == 'tissue_ontology_term_id_0'
                        and coord[3] == 'organism_ontology_term_id_0'):
                    attr_vals['nnz'][i] = 1
                    attr_vals['n_cells'][i] = 10
                    attr_vals['sum'][i] = 3.0
            return attr_vals

        with create_temp_cube(dim_size=2, attr_vals_fn=attr_values) as cube:
            qry = cube.query()

            res = qry.df['feature_id_0', 'tissue_ontology_term_id_0', 'organism_ontology_term_id_0']
            assert all(res['nnz'] == 1)
            assert all(res['n_cells'] == 10)
            assert all(res['sum'] == 3.0)

            res = qry.df['feature_id_1', 'tissue_ontology_term_id_0', 'organism_ontology_term_id_0']
            assert all(res['nnz'] == 0)
            assert all(res['n_cells'] == 0)
            assert all(res['sum'] == 0)
