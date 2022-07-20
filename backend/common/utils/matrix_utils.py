import logging

import numpy as np


def is_matrix_sparse(matrix: np.ndarray, sparse_threshold):
    """
    Returns whether `matrix` is sparse or not (i.e. dense). This is determined by figuring out whether the matrix has
    a sparsity percentage below the sparse_threshold, returning the number of non-zeros encountered and number of
    elements evaluated.  This function may return before evaluating the whole matrix if it can be determined that matrix
    is not sparse enough.
    """

    if sparse_threshold == 100.0:
        return True
    if sparse_threshold == 0.0:
        return False

    total_number_of_rows = matrix.shape[0]
    total_number_of_columns = matrix.shape[1]
    total_number_of_matrix_elements = total_number_of_rows * total_number_of_columns

    # For efficiency, we count the number of non-zero elements in chunks of the matrix at a time until we hit the
    # maximum number of non zero values allowed before the matrix is deemed "dense." This allows the function the
    # quit early for large dense matrices.
    row_stride = min(int(np.power(10, np.around(np.log10(1e9 / total_number_of_columns)))), 10_000)

    maximum_number_of_non_zero_elements_in_matrix = int(
        total_number_of_rows * total_number_of_columns * sparse_threshold / 100
    )
    number_of_non_zero_elements = 0

    for start_row_index in range(0, total_number_of_rows, row_stride):
        end_row_index = min(start_row_index + row_stride, total_number_of_rows)

        matrix_subset = matrix[start_row_index:end_row_index, :]
        if not isinstance(matrix_subset, np.ndarray):
            matrix_subset = matrix_subset.toarray()

        number_of_non_zero_elements += np.count_nonzero(matrix_subset)
        if number_of_non_zero_elements > maximum_number_of_non_zero_elements_in_matrix:
            if end_row_index != total_number_of_rows:
                percentage_of_non_zero_elements = (
                    100 * number_of_non_zero_elements / (end_row_index * total_number_of_columns)
                )
                logging.info(
                    f"Matrix is not sparse. Percentage of non-zero elements (estimate): "
                    f"{percentage_of_non_zero_elements:6.2f}"
                )
            else:
                percentage_of_non_zero_elements = 100 * number_of_non_zero_elements / total_number_of_matrix_elements
                logging.info(
                    f"Matrix is not sparse. Percentage of non-zero elements (exact): "
                    f"{percentage_of_non_zero_elements:6.2f}"
                )
            return False

    is_sparse = (100.0 * number_of_non_zero_elements / total_number_of_matrix_elements) < sparse_threshold
    return is_sparse
