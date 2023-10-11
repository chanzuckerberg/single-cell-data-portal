import tempfile
import unittest

import numpy as np

from server.common.fbs.matrix import decode_matrix_fbs, encode_matrix_fbs
from server.compute import diffexp_cxg
from server.compute.diffexp_cxg import diffexp_ttest
from server.dataset.matrix_loader import DataLoader
from server.tests import FIXTURES_ROOT
from server.tests.unit import app_config


class DiffExpTest(unittest.TestCase):
    """Tests the diffexp returns the expected results for one test case, using different
    adaptor types and different algorithms."""

    def load_dataset(self, path, extra_server_config=None, extra_dataset_config=None):
        extra_server_config = {} if extra_server_config is None else extra_server_config
        extra_dataset_config = {} if extra_dataset_config is None else extra_dataset_config
        extra_dataset_config["X_approximate_distribution"] = "normal"  # hardwired for now
        config = app_config(extra_server_config=extra_server_config, extra_dataset_config=extra_dataset_config)
        loader = DataLoader(location=path, app_config=config)
        adaptor = loader.open()
        return adaptor

    def get_mask(self, adaptor, start, stride):
        """Simple function to return a mask or rows"""
        rows = adaptor.get_shape()[0]
        sel = list(range(start, rows, stride))
        mask = np.zeros(rows, dtype=bool)
        mask[sel] = True
        return mask

    def compare_diffexp_results(self, results, expects):
        self.assertEqual(len(results), len(expects))
        for result, expect in zip(results, expects):
            self.assertEqual(result[0], expect[0])
            self.assertTrue(np.isclose(result[1], expect[1], 1e-6, 1e-4))
            self.assertTrue(np.isclose(result[2], expect[2], 1e-6, 1e-4))
            self.assertTrue(np.isclose(result[3], expect[3], 1e-6, 1e-4))

    def check_1_10_2_10(self, results):
        """Checks the results for a specific set of rows selections"""

        positive_expects = [
            [1712, 0.24104056, 0.0051788902660723345, 1.0],
            [1575, 0.2615018, 0.007830310753043345, 1.0],
            [693, 0.23106655, 0.008715846769131548, 1.0],
            [916, 0.2395215, 0.009080596532247588, 1.0],
            [77, 0.22927025, 0.010070392939027756, 1.0],
            [782, 0.20581803, 0.010161745218916036, 1.0],
            [913, 0.23841085, 0.010782030711612685, 1.0],
            [910, 0.21493295, 0.014596411069229197, 1.0],
            [1727, 0.21911663, 0.015168372104237176, 1.0],
            [1443, 0.19814226, 0.015337080567465522, 1.0],
        ]
        negative_expects = [
            [956, -0.29662406, 0.0008649321884808977, 1.0],
            [1124, -0.2607333, 0.0011717216548271284, 1.0],
            [1809, -0.24854594, 0.0019304405196777848, 1.0],
            [1754, -0.24683577, 0.005691734062127954, 1.0],
            [948, -0.18708363, 0.006622111055981219, 1.0],
            [1810, -0.2172082, 0.007055917428377063, 1.0],
            [779, -0.21150622, 0.007202934422407284, 1.0],
            [576, -0.19008157, 0.008272092578813124, 1.0],
            [538, -0.21803819, 0.01062259019889307, 1.0],
            [436, -0.2100364, 0.01127515110543434, 1.0],
        ]

        self.compare_diffexp_results(results["positive"], positive_expects)
        self.compare_diffexp_results(results["negative"], negative_expects)

    def get_X_col(self, adaptor, cols):
        varmask = np.zeros(adaptor.get_shape()[1], dtype=bool)
        varmask[cols] = True
        return adaptor.get_X_array(None, varmask)

    def test_cxg_default(self):
        """Test a cxg adaptor with its default diffexp algorithm (diffexp_cxg)"""
        for dataset in ["pbmc3k.cxg", "pbmc3k_sparse.cxg"]:
            with self.subTest(dataset=dataset):
                adaptor = self.load_dataset(f"{FIXTURES_ROOT}/{dataset}")
                maskA = self.get_mask(adaptor, 1, 10)
                maskB = self.get_mask(adaptor, 2, 10)

                # run it through the adaptor
                results = adaptor.compute_diffexp_ttest(maskA, maskB, 10)
                self.check_1_10_2_10(results)

                # check list (not mask) variant
                results = adaptor.compute_diffexp_ttest(maskA.nonzero()[0], maskB.nonzero()[0], 10, selector_lists=True)
                self.check_1_10_2_10(results)

                # run it directly

                results = diffexp_ttest(adaptor, maskA, maskB, 10)
                self.check_1_10_2_10(results)

                results = diffexp_ttest(adaptor, maskA.nonzero()[0], maskB.nonzero()[0], 10, selector_lists=True)
                self.check_1_10_2_10(results)

    def test_cxg_sparse(self):
        adaptor_sparse = self.load_dataset(
            f"{FIXTURES_ROOT}/diffexp/sparse_no_col_shift.cxg",
        )
        adaptor_dense = self.load_dataset(f"{FIXTURES_ROOT}/diffexp/dense_no_col_shift.cxg")

        self.sparse_diffexp(adaptor_dense, adaptor_sparse)

    def sparse_diffexp(self, adaptor_dense, adaptor_sparse):
        with tempfile.TemporaryDirectory():
            maskA = self.get_mask(adaptor_dense, 1, 10)
            maskB = self.get_mask(adaptor_dense, 2, 10)

            diffexp_results_sparse = diffexp_cxg.diffexp_ttest(adaptor_sparse, maskA, maskB, 10)
            diffexp_results_dense = diffexp_cxg.diffexp_ttest(adaptor_dense, maskA, maskB, 10)

            self.compare_diffexp_results(diffexp_results_dense["positive"], diffexp_results_sparse["positive"])
            self.compare_diffexp_results(diffexp_results_dense["negative"], diffexp_results_sparse["negative"])

            topcols_pos = np.array([x[0] for x in diffexp_results_dense["positive"]])
            topcols_neg = np.array([x[0] for x in diffexp_results_dense["negative"]])
            topcols = np.concatenate((topcols_pos, topcols_neg))

            cols_sparse = self.get_X_col(adaptor_sparse, topcols)
            cols_dense = self.get_X_col(adaptor_dense, topcols)

            assert cols_dense.shape[0] == adaptor_sparse.get_shape()[0]
            assert cols_dense.shape[1] == len(diffexp_results_dense["positive"]) + len(
                diffexp_results_dense["negative"]
            )

            def convert(mat, cols):
                return decode_matrix_fbs(encode_matrix_fbs(mat, col_idx=cols)).to_numpy()

            cols_sparse = convert(cols_sparse, topcols)
            cols_dense = convert(cols_dense, topcols)

            x = adaptor_sparse.get_X_array()
            assert x.shape == adaptor_sparse.get_shape()

            for row in range(cols_dense.shape[0]):
                for col in range(cols_dense.shape[1]):
                    vsparse = cols_sparse[row][col]
                    vdense = cols_dense[row][col]
                    self.assertTrue(np.isclose(vdense, vsparse, 1e-6, 1e-6))
