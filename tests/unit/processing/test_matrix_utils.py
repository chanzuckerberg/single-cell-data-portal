import logging
from unittest.mock import Mock

import numpy as np
import pytest
from anndata import AnnData
from scipy.sparse import coo_matrix

from backend.layers.processing.utils.matrix_utils import enforce_canonical_format, is_matrix_sparse

LOGGER = logging.getLogger("matrix_utils")
LOGGER.propagate = True


class TestMatrixUtils:
    def test__is_matrix_sparse__zero_and_one_hundred_percent_threshold(self):
        matrix = np.array([1, 2, 3])

        assert not is_matrix_sparse(matrix, 0)
        assert is_matrix_sparse(matrix, 100)

    def test__is_matrix_sparse__partially_populated_sparse_matrix_returns_true(self):
        matrix = np.zeros([3, 4])
        matrix[2][3] = 1.0
        matrix[1][1] = 2.2

        assert is_matrix_sparse(matrix, 50)

    def test__is_matrix_sparse__partially_populated_dense_matrix_returns_false(self):
        matrix = np.zeros([2, 2])
        matrix[0][0] = 1.0
        matrix[0][1] = 2.2
        matrix[1][1] = 3.7

        assert not is_matrix_sparse(matrix, 50)

    def test__is_matrix_sparse__giant_matrix_returns_false_early(self, caplog):
        caplog.set_level(logging.INFO)
        matrix = np.ones([20000, 20])

        assert not is_matrix_sparse(matrix, 1)

        # Because the function returns early a log will output the _estimate_ instead of the _exact_ percentage of
        # non-zero elements in the matrix.
        assert "Percentage of non-zero elements (estimate)" in caplog.text

    def test__is_matrix_sparse_with_column_shift_encoding__giant_matrix_returns_false_early(self, caplog):
        caplog.set_level(logging.INFO)
        matrix = np.random.rand(20000, 20)

        assert not is_matrix_sparse(matrix, 1)

        # Because the function returns early a log will output the _estimate_ instead of the _exact_ percentage of
        # non-zero elements in the matrix.
        assert "Percentage of non-zero elements (estimate)" in caplog.text


@pytest.fixture
def noncanonical_matrix():
    array = np.array([[1, 0, 1], [3, 2, 3], [4, 5, 4]])
    return coo_matrix((array[0], (array[1], array[2])))


@pytest.fixture
def canonical_adata():
    return Mock(X=Mock(has_canonical_format=True))


class TestEnforceCanonical:
    def test_adata_with_noncanonical_X_and_raw_X(self, noncanonical_matrix, caplog):
        assert noncanonical_matrix.has_canonical_format is False
        adata = AnnData(noncanonical_matrix)
        enforce_canonical_format(adata)
        assert adata.X.has_canonical_format is True
        assert "noncanonical data found in X; converting to canonical format using sum_duplicates." in caplog.text

    def test_adata_with_noncanonical_raw_X(self, noncanonical_matrix, caplog):
        caplog.set_level(logging.WARNING)
        assert noncanonical_matrix.has_canonical_format is False
        adata = AnnData(raw=AnnData(noncanonical_matrix))
        enforce_canonical_format(adata.raw)
        assert adata.raw.X.has_canonical_format is True
        assert "noncanonical data found in X; converting to canonical format using sum_duplicates." in caplog.text

    def test_adata_with_canonical_X(self, canonical_adata, caplog):
        caplog.set_level(logging.WARNING)
        enforce_canonical_format(canonical_adata)
        assert canonical_adata.X.has_canonical_format is True
        assert "noncanonical data found in X; converting to canonical format using sum_duplicates." not in caplog.text
