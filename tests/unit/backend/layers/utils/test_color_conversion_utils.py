import unittest

import anndata

from backend.common.utils.color_conversion_utils import (
    ColorFormatException,
    convert_anndata_category_colors_to_cxg_category_colors,
    convert_color_to_hex_format,
)
from tests.unit.backend.fixtures.environment_setup import fixture_file_path


class TestColorConversionUtils(unittest.TestCase):
    """Test color conversion helper functions"""

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

    def test_term_vs_label_colors(self):
        expected_label_colors = {
            "assay": {"10x 3' v2": "#ffb6c1", "10x 3' v3": "#ffc0cb"},
            "cell_type": {"mature B cell": "#ffb6c1", "plasma cell": "#ffc0cb"},
            "development_stage": {"human adult stage": "#ff69b4", "mature stage": "#ffc0cb"},
            "donor_id": {"C41": "#ff0000", "C58": "#ffa500", "C70": "#ffff00", "C72": "#008000"},
            "sex": {"female": "#ffc0cb", "male": "#0000ff"},
            "suspension_type": {"cell": "#0000ff", "nucleus": "#800080"},
        }
        adata_with_term_colors = anndata.read_h5ad(fixture_file_path("liver_with_term_colors.h5ad"))
        self.assertDictEqual(
            convert_anndata_category_colors_to_cxg_category_colors(adata_with_term_colors), expected_label_colors
        )

        adata_with_label_colors = anndata.read_h5ad(fixture_file_path("liver_with_label_colors.h5ad"))
        self.assertDictEqual(
            convert_anndata_category_colors_to_cxg_category_colors(adata_with_label_colors), expected_label_colors
        )

    def _get_h5ad(self):
        return anndata.read_h5ad(fixture_file_path("pbmc3k.h5ad"))
