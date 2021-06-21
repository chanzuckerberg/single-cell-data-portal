import unittest

import anndata

from backend.corpora.common.utils.color_conversion_utils import (
    convert_color_to_hex_format,
    convert_anndata_category_colors_to_cxg_category_colors,
)
from backend.corpora.common.utils.exceptions import ColorFormatException
from tests.unit.backend import PROJECT_ROOT


class ColorsTest(unittest.TestCase):
    """ Test color helper functions """

    def test_convert_color_to_hex_format(self):
        self.assertEqual(convert_color_to_hex_format("wheat"), "#f5deb3")
        self.assertEqual(convert_color_to_hex_format("WHEAT"), "#f5deb3")
        self.assertEqual(convert_color_to_hex_format((245, 222, 179)), "#f5deb3")
        self.assertEqual(convert_color_to_hex_format([245, 222, 179]), "#f5deb3")
        self.assertEqual(convert_color_to_hex_format("#f5deb3"), "#f5deb3")
        self.assertEqual(
            convert_color_to_hex_format([0.9607843137254902, 0.8705882352941177, 0.7019607843137254]), "#f5deb3"
        )
        for bad_input in ["foo", "BAR", "#AABB", "#AABBCCDD", "#AABBGG", (1, 2), [1, 2], (1, 2, 3, 4), [1, 2, 3, 4]]:
            with self.assertRaises(ColorFormatException):
                convert_color_to_hex_format(bad_input)

    def test_anndata_colors_to_cxg_colors(self):
        # test standard behavior
        adata = self._get_h5ad()
        expected_pbmc3k_colors = {
            "louvain": {
                "B cells": "#2ca02c",
                "CD14+ Monocytes": "#ff7f0e",
                "CD4 T cells": "#1f77b4",
                "CD8 T cells": "#d62728",
                "Dendritic cells": "#e377c2",
                "FCGR3A+ Monocytes": "#8c564b",
                "Megakaryocytes": "#bcbd22",
                "NK cells": "#9467bd",
            }
        }
        self.assertEqual(convert_anndata_category_colors_to_cxg_category_colors(adata), expected_pbmc3k_colors)

        # test that invalid color formats raise an exception
        adata.uns["louvain_colors"][0] = "#NOTCOOL"
        with self.assertRaises(ColorFormatException):
            convert_anndata_category_colors_to_cxg_category_colors(adata)

        # test that colors without a matching obs category are skipped
        adata = self._get_h5ad()
        del adata.obs["louvain"]
        self.assertEqual(convert_anndata_category_colors_to_cxg_category_colors(adata), {})

    def _get_h5ad(self):
        return anndata.read_h5ad(f"{PROJECT_ROOT}/tests/unit/backend/corpora/fixtures/pbmc3k.h5ad")
