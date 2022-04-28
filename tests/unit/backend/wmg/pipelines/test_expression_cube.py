import tempfile
import unittest

import tiledb

from backend.wmg.data.schemas.cube_schema import expression_summary_schema
from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_CUBE_NAME
from backend.wmg.data.tiledb import create_ctx


class TestCellExpressionSummaryCubeETL(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        super().setUp(cls)
        cls.tmp_dir = tempfile.mkdtemp()
        expected_expression_summary_uri = f"{cls.tmp_dir}/fixtures/{EXPRESSION_SUMMARY_CUBE_NAME}"
        # create empty expression summary cube
        ctx = create_ctx()

        with tiledb.scope_ctx(ctx):
            tiledb.Array.create(expected_expression_summary_uri, expression_summary_schema, overwrite=True)
