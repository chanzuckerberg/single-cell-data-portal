import json
import pickle
import unittest
from os import mkdir, path
from shutil import rmtree
from uuid import uuid4

import numpy as np
import tiledb
from pandas import Categorical, DataFrame, Series

from backend.common.utils.cxg_generation_utils import (
    convert_dataframe_to_cxg_array,
    convert_dictionary_to_cxg_group,
    convert_matrices_to_cxg_arrays,
    convert_ndarray_to_cxg_dense_array,
    convert_uns_to_cxg_group,
)
from tests.unit.backend.fixtures.environment_setup import fixture_file_path


class TestCxgGenerationUtils(unittest.TestCase):
    def setUp(self):
        self.testing_cxg_temp_directory = fixture_file_path(str(uuid4()))
        mkdir(self.testing_cxg_temp_directory)

    def tearDown(self):
        if path.isdir(self.testing_cxg_temp_directory):
            rmtree(self.testing_cxg_temp_directory)

    def test__convert_dictionary_to_cxg_group__writes_successfully(self):
        random_dictionary = {"cookies": "chocolate_chip", "brownies": "chocolate", "cake": "double chocolate"}
        dictionary_name = "spatial"
        expected_array_directory = f"{self.testing_cxg_temp_directory}/{dictionary_name}"

        convert_dictionary_to_cxg_group(
            self.testing_cxg_temp_directory, random_dictionary, group_metadata_name=dictionary_name
        )

        array = tiledb.open(expected_array_directory)
        actual_stored_metadata = dict(array.meta.items())

        self.assertTrue(path.isdir(expected_array_directory))
        self.assertTrue(isinstance(array, tiledb.DenseArray))
        self.assertEqual(random_dictionary, actual_stored_metadata)

    def test__convert_uns_to_cxg_group__writes_successfully(self):
        random_dictionary = {
            "spatial": {
                "abcd": {
                    "images": {
                        "hires": "123",
                        "fullres": [],
                    },
                    "scalefactors": {
                        "spot_diameter_fullres": "123",
                        "tissue_hires_scalef": "123",
                    },
                }
            }
        }
        dictionary_name = "uns"
        expected_array_directory = f"{self.testing_cxg_temp_directory}/{dictionary_name}"
        convert_uns_to_cxg_group(
            self.testing_cxg_temp_directory, random_dictionary, group_metadata_name=dictionary_name
        )
        array = tiledb.open(expected_array_directory)
        actual_stored_metadata = dict(array.meta.items())

        self.assertTrue(path.isdir(expected_array_directory))
        self.assertTrue(isinstance(array, tiledb.DenseArray))
        self.assertEqual(random_dictionary["spatial"], pickle.loads(actual_stored_metadata["spatial"]))

    def test__convert_dataframe_to_cxg_array__writes_successfully(self):
        random_int_category = Series(data=[3, 1, 2, 4], dtype=np.int64)
        random_bool_category = Series(data=[True, True, False, True], dtype=np.bool_)
        random_nan_category = Categorical(Series(data=["A", "B", np.nan, "D"]))
        random_dataframe_name = f"random_dataframe_{uuid4()}"
        random_dataframe = DataFrame(
            data={
                "int_category": random_int_category,
                "bool_category": random_bool_category,
                "nan_category": random_nan_category,
            }
        )

        convert_dataframe_to_cxg_array(
            self.testing_cxg_temp_directory, random_dataframe_name, random_dataframe, "int_category", tiledb.Ctx()
        )

        expected_array_directory = f"{self.testing_cxg_temp_directory}/{random_dataframe_name}"
        expected_array_metadata = {
            "cxg_schema": json.dumps(
                {
                    "int_category": {"type": "int32"},
                    "bool_category": {"type": "boolean"},
                    "nan_category": {"type": "categorical", "categories": ["A", "B", "D", "nan"]},
                    "index": "int_category",
                }
            )
        }

        actual_stored_dataframe_array = tiledb.open(expected_array_directory)
        actual_stored_dataframe_metadata = dict(actual_stored_dataframe_array.meta.items())

        self.assertTrue(path.isdir(expected_array_directory))
        self.assertTrue(isinstance(actual_stored_dataframe_array, tiledb.DenseArray))
        self.assertDictEqual(expected_array_metadata, actual_stored_dataframe_metadata)
        self.assertTrue((actual_stored_dataframe_array[0:4]["int_category"] == random_int_category.to_numpy()).all())
        self.assertTrue((actual_stored_dataframe_array[0:4]["bool_category"] == random_bool_category.to_numpy()).all())

    def test__convert_ndarray_to_cxg_dense_array__writes_successfully(self):
        ndarray = np.random.rand(3, 2)
        ndarray_name = f"{self.testing_cxg_temp_directory}/awesome_ndarray_{uuid4()}"

        convert_ndarray_to_cxg_dense_array(ndarray_name, ndarray, tiledb.Ctx())

        actual_stored_array = tiledb.open(ndarray_name)

        self.assertTrue(path.isdir(ndarray_name))
        self.assertTrue(isinstance(actual_stored_array, tiledb.DenseArray))
        self.assertTrue((actual_stored_array[:, :] == ndarray).all())

    def test__convert_matrices_to_cxg_arrays__dense_array_writes_successfully(self):
        matrix = np.float32(np.random.rand(3, 2))
        matrix_name = f"{self.testing_cxg_temp_directory}/awesome_matrix_{uuid4()}"

        convert_matrices_to_cxg_arrays(matrix_name, matrix, False, tiledb.Ctx())
        actual_stored_array = tiledb.open(matrix_name)
        self.assertTrue(path.isdir(matrix_name))
        self.assertTrue(isinstance(actual_stored_array, tiledb.DenseArray))
        self.assertTrue((actual_stored_array[:, :] == matrix).all())

    def test__convert_matrices_to_cxg_arrays__sparse_array_only_store_nonzeros_empty_array(self):
        matrix = np.zeros([3, 2])
        matrix_name = f"{self.testing_cxg_temp_directory}/awesome_zero_matrix_{uuid4()}"

        convert_matrices_to_cxg_arrays(matrix_name, matrix, True, tiledb.Ctx())

        for suffix in ["r", "c"]:
            actual_stored_array = tiledb.open(matrix_name + suffix)
            self.assertTrue(path.isdir(matrix_name + suffix))
            self.assertTrue(isinstance(actual_stored_array, tiledb.SparseArray))
            self.assertTrue(actual_stored_array[:][""].size == 0)

    def test__convert_matrices_to_cxg_arrays__sparse_array_only_store_nonzeros(self):
        matrix = np.zeros([3, 3])
        matrix[0, 0] = 1
        matrix[1, 1] = 1
        matrix[2, 2] = 2
        matrix_name = f"{self.testing_cxg_temp_directory}/awesome_sparse_matrix_{uuid4()}"

        convert_matrices_to_cxg_arrays(matrix_name, matrix, True, tiledb.Ctx())

        def get_value_at_coord(array, coord, attr):
            x, y = coord
            return array[x][""][array[x][attr] == y][0]

        for suffix, attr_dim in zip(["r", "c"], ["var", "obs"]):
            actual_stored_array = tiledb.open(matrix_name + suffix)
            self.assertTrue(path.isdir(matrix_name + suffix))
            self.assertTrue(isinstance(actual_stored_array, tiledb.SparseArray))
            self.assertTrue(get_value_at_coord(actual_stored_array, (0, 0), attr_dim) == 1)
            self.assertTrue(get_value_at_coord(actual_stored_array, (1, 1), attr_dim) == 1)
            self.assertTrue(get_value_at_coord(actual_stored_array, (2, 2), attr_dim) == 2)
            self.assertTrue(actual_stored_array[:][""].size == 3)
