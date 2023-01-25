import unittest
import numpy as np
from backend.wmg.api.rollup import rollup_across_cell_type_descendants


class RollupExpressionsAcrossCellTypesTest(unittest.TestCase):
    # test that the rollup function works as expected for input 2d arrays
    def test__expression_rollup_2d_arrays(self):
        # second cell type is descendant of first cell type
        # fourth cell type is descendant of third cell type
        cell_types = ["CL:0000786", "CL:0000986", "CL:0000980", "CL:0001202"]
        np.random.seed(0)
        exprs = [np.random.rand(len(cell_types), 10) for _ in range(3)]
        exprs_rollup = rollup_across_cell_type_descendants(cell_types, exprs)
        expected = []
        for array in exprs:
            expected_array = array.copy()
            expected_array[0] = expected_array[0] + expected_array[1]
            expected_array[2] = expected_array[2] + expected_array[3]
            expected.append(expected_array)
        assert np.all([np.array_equal(exprs_rollup[i], expected[i]) for i in range(len(expected))])

    # test that the rollup function works as expected for input 1d arrays
    def test__expression_rollup_1d_arrays(self):
        # second cell type is descendant of first cell type
        # fourth cell type is descendant of third cell type
        cell_types = ["CL:0000786", "CL:0000986", "CL:0000980", "CL:0001202"]

        np.random.seed(0)
        exprs = [np.random.rand(len(cell_types)) for _ in range(3)]
        exprs_rollup = rollup_across_cell_type_descendants(cell_types, exprs)
        expected = []
        for array in exprs:
            expected_array = array.copy()
            expected_array[0] = expected_array[0] + expected_array[1]
            expected_array[2] = expected_array[2] + expected_array[3]
            expected.append(expected_array)
        assert np.all([np.array_equal(exprs_rollup[i], expected[i]) for i in range(len(expected))])

    def test__expression_rollup_no_descendants_overlap(self):
        # these cell types are not descendants of each other
        cell_types = ["CL:0000786", "CL:0000980"]
        np.random.seed(0)
        exprs = [np.random.rand(len(cell_types), 10) for _ in range(3)]
        exprs_rollup = rollup_across_cell_type_descendants(cell_types, exprs)
        expected = []
        for array in exprs:
            expected.append(array)
        assert np.all([np.array_equal(exprs_rollup[i], expected[i]) for i in range(len(expected))])
