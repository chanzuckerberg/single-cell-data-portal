"""This module tests the low level functions encapsulating the roll up operation across the cell type ontology.

In detail, this module tests the public and private functions defined in `backend.common.utils.rollup` module.
"""

import unittest

import numpy as np
import pandas as pd

from backend.common.utils.rollup import (
    are_cell_types_colinear,
    are_cell_types_not_redundant_nodes,
    rollup_across_cell_type_descendants,
)


class TestLowLevelRollupFunctionsTraversingCellTypeLineage(unittest.TestCase):
    # test that the rollup function works as expected
    def test__expression_rollup_across_cell_type_descendants(self):
        # second cell type is descendant of first cell type
        # fourth cell type is descendant of third cell type
        cell_types = ["CL:0000786", "CL:0000986", "CL:0000980", "CL:0001202"]
        df = pd.DataFrame()
        df["cell_type_ontology_term_id"] = cell_types
        np.random.seed(0)
        exprs = np.random.rand(len(cell_types), 10)
        for i in range(exprs.shape[1]):
            df[i] = exprs[:, i]

        df_rollup = rollup_across_cell_type_descendants(df)
        df_expected = df.copy()
        expected_exprs = exprs.copy()
        expected_exprs[0] = exprs[0] + exprs[1]
        expected_exprs[2] = exprs[2] + exprs[3]

        for i in range(expected_exprs.shape[1]):
            df_expected[i] = expected_exprs[:, i]
        assert np.all(df_expected == df_rollup)

    def test__expression_rollup_no_descendants_overlap(self):
        # these cell types are not descendants of each other
        cell_types = ["CL:0000786", "CL:0000980"]
        df = pd.DataFrame()
        df["cell_type_ontology_term_id"] = cell_types
        np.random.seed(0)
        exprs = np.random.rand(len(cell_types), 10)
        for i in range(exprs.shape[1]):
            df[i] = exprs[:, i]

        df_rollup = rollup_across_cell_type_descendants(df)
        assert np.all(df == df_rollup)

    def test__cell_types_in_same_lineage_are_colinear(self):
        # first and second pairs of cell types are in the same lineage
        # third pair of cell types are not in the same lineage
        cell_type_pairs = [["CL:0000786", "CL:0000986"], ["CL:0000980", "CL:0001202"], ["CL:0000786", "CL:0000980"]]
        expected_colinearity = [True, True, False]
        for cell_types, expected in zip(cell_type_pairs, expected_colinearity):
            a, b = cell_types
            assert are_cell_types_colinear(a, b) == expected

    def test__cell_types_are_not_redundant(self):
        cell_counts = {
            "CL:0000127;;UBERON:0000955--HANCESTRO:0005": 10,
            "CL:0000127;;UBERON:0000955--HANCESTRO:0006": 10,
            "CL:0000127;;UBERON:0000955--HANCESTRO:0008": 50,
            "CL:0000127;;UBERON:0000955--multiethnic": 70,
            "CL:0000127;;UBERON:0000955--unknown": 410,
            "CL:0000644;;UBERON:0000955--unknown": 70,
            "CL:0002605;;UBERON:0000955--HANCESTRO:0005": 10,
            "CL:0002605;;UBERON:0000955--HANCESTRO:0008": 30,
            "CL:0002605;;UBERON:0000955--multiethnic": 40,
            "CL:0002627;;UBERON:0000955--HANCESTRO:0006": 10,
            "CL:0002627;;UBERON:0000955--HANCESTRO:0008": 20,
            "CL:0002627;;UBERON:0000955--multiethnic": 30,
            "CL:0002627;;UBERON:0000955--unknown": 40,
        }
        cell_type_groups = np.array(
            [
                "CL:0000127;;UBERON:0000955--HANCESTRO:0005",
                "CL:0000127;;UBERON:0000955--HANCESTRO:0006",
                "CL:0000127;;UBERON:0000955--HANCESTRO:0008",
                "CL:0000127;;UBERON:0000955--multiethnic",
                "CL:0000127;;UBERON:0000955--unknown",
                "CL:0000644;;UBERON:0000955--unknown",
                "CL:0002605;;UBERON:0000955--HANCESTRO:0005",
                "CL:0002605;;UBERON:0000955--HANCESTRO:0008",
                "CL:0002605;;UBERON:0000955--multiethnic",
                "CL:0002627;;UBERON:0000955--HANCESTRO:0006",
                "CL:0002627;;UBERON:0000955--HANCESTRO:0008",
                "CL:0002627;;UBERON:0000955--multiethnic",
                "CL:0002627;;UBERON:0000955--unknown",
            ],
            dtype="<U46",
        )

        result = are_cell_types_not_redundant_nodes(cell_type_groups, cell_counts)
        expected = [False, False, True, True, True, True, True, True, True, True, True, True, True]
        self.assertEqual(result, expected)

    def test__cell_types_are_not_redundant_linear_chain(self):
        cell_counts = {
            "CL:0000127": 10,
            "CL:0000683": 50,
            "CL:0002085": 50,
        }
        cell_type_groups = np.array(["CL:0000127", "CL:0000683", "CL:0002085"])

        result = are_cell_types_not_redundant_nodes(cell_type_groups, cell_counts)
        expected = [True, False, True]
        self.assertEqual(result, expected)
