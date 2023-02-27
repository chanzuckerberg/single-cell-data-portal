import unittest

import numpy as np

from backend.wmg.data.schemas.cube_schema_default import (
    expression_summary_logical_dims as expression_summary_default_logical_dims,
)
from tests.unit.backend.wmg.fixtures.test_snapshot import load_realistic_test_snapshot

TEST_SNAPSHOT = "realistic-test-snapshot"

# creates the expression summary default cube by group-by and summing the expression summary cube
# and checks to makes sure that the result is the same as the existing fixture (expression_summary_default)


class ExpressionSummaryDefaultCubeCreation(unittest.TestCase):
    def test__expression_summary_and_expression_summary_default_cube_equality(self):
        with load_realistic_test_snapshot(TEST_SNAPSHOT) as snapshot:
            expression_summary_df = snapshot.expression_summary_cube.df[:]
            expression_summary_default_df = snapshot.expression_summary_default_cube.df[:]

            # generate expression summary default cube by group-by and summing the expression summary cube
            test_expression_summary_default_df = (
                expression_summary_df.groupby(expression_summary_default_logical_dims).sum().reset_index()
            )

            # check that the two cubes are equal
            df1 = expression_summary_default_df.sort_values(expression_summary_default_logical_dims)
            df2 = test_expression_summary_default_df.sort_values(expression_summary_default_logical_dims)
            rec1 = df1.to_dict(orient="records")
            rec2 = df2.to_dict(orient="records")
            # return expression_summary_default_df,test_expression_summary_default_df,rec1,rec2
            [np.testing.assert_equal(rec1[i], rec2[i]) for i in range(len(rec1))]
