import json
import unittest
from os import path, remove
from shutil import rmtree
from uuid import uuid4

import anndata
import numpy as np
import tiledb
from pandas import Categorical, DataFrame, Series

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.processing.h5ad_data_file import H5ADDataFile
from tests.unit.backend.fixtures.environment_setup import fixture_file_path


class TestH5ADDataFile(unittest.TestCase):
    def setUp(self):
        self.sample_anndata = self._create_sample_anndata_dataset()
        self.sample_h5ad_filename = self._write_anndata_to_file(self.sample_anndata)

        self.sample_output_directory = path.splitext(self.sample_h5ad_filename)[0] + ".cxg"

        def mock_config_fn(name):
            if name == "schema_4_feature_flag":
                return "True"

        self.mock_config = unittest.mock.patch(
            "backend.common.corpora_config.CorporaConfig.__getattr__", side_effect=mock_config_fn
        )
        self.mock_config.start()

    def tearDown(self):
        if self.sample_h5ad_filename:
            remove(self.sample_h5ad_filename)

        if path.isdir(self.sample_output_directory):
            rmtree(self.sample_output_directory)

        self.mock_config.stop()

    def test__create_h5ad_data_file__non_h5ad_raises_exception(self):
        non_h5ad_filename = "my_fancy_dataset.csv"

        with self.assertRaises(Exception) as exception_context:
            H5ADDataFile(non_h5ad_filename)

        self.assertIn("File must be an H5AD", str(exception_context.exception))

    def test__create_h5ad_data_file__reads_anndata_successfully(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)

        self.assertTrue((h5ad_file.anndata.X == self.sample_anndata.X).all())
        self.assertEqual(
            h5ad_file.anndata.obs.sort_index(inplace=True), self.sample_anndata.obs.sort_index(inplace=True)
        )
        self.assertEqual(
            h5ad_file.anndata.var.sort_index(inplace=True), self.sample_anndata.var.sort_index(inplace=True)
        )

        for key in h5ad_file.anndata.obsm:
            self.assertIn(key, self.sample_anndata.obsm.keys())
            self.assertTrue((h5ad_file.anndata.obsm[key] == self.sample_anndata.obsm[key]).all())

        for key in self.sample_anndata.obsm:
            self.assertIn(key, h5ad_file.anndata.obsm.keys())
            self.assertTrue((h5ad_file.anndata.obsm[key] == self.sample_anndata.obsm[key]).all())

    def test__create_h5ad_data_file__copies_index_of_obs_and_var_to_column(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)

        # The automatic name chosen for the index should be "name_0"
        self.assertNotIn("name_0", self.sample_anndata.obs.columns)
        self.assertIn("name_0", h5ad_file.obs.columns)

        self.assertNotIn("name_0", self.sample_anndata.var.columns)
        self.assertIn("name_0", h5ad_file.var.columns)

    def test__create_h5ad_data_file__no_copy_if_var_index_name_already_in_columns(self):
        h5ad_file = H5ADDataFile(
            self.sample_h5ad_filename,
            var_index_column_name="int_category",
        )
        self.assertNotIn("name_0", h5ad_file.var.columns)
        self.assertIn("int_category", h5ad_file.var.columns)

    def test__create_h5ad_data_file__copy_if_var_index_name_not_in_columns(self):
        h5ad_file = H5ADDataFile(
            self.sample_h5ad_filename,
            var_index_column_name="int_category",
        )
        # set index that is not already in columns
        anndata_object = anndata.read_h5ad(self.sample_h5ad_filename)
        anndata_object.var.reset_index(drop=True, inplace=True)

        var_with_index_col = h5ad_file.transform_dataframe_index_into_column(anndata_object.var, "var", "int_category")

        self.assertIn("name_0", var_with_index_col.columns)
        self.assertIn("int_category", var_with_index_col.columns)

    def test__create_h5ad_data_file__var_index_name_specified_not_unique_raises_exception(self):

        with self.assertRaises(Exception) as exception_context:
            H5ADDataFile(
                self.sample_h5ad_filename,
                var_index_column_name="bool_category",
            )

        self.assertIn("Please prepare data to contain unique values", str(exception_context.exception))

    def test__create_h5ad_data_file__obs_and_var_index_names_specified_doesnt_exist_raises_exception(self):
        with self.assertRaises(Exception) as exception_context:
            H5ADDataFile(
                self.sample_h5ad_filename,
                var_index_column_name="i_dont_exist",
            )

        self.assertIn("does not exist", str(exception_context.exception))

    def test__to_cxg__simple_anndata_no_corpora_and_sparse(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 100)

        self._validate_cxg_and_h5ad_content_match(self.sample_h5ad_filename, self.sample_output_directory, True)

    def test__to_cxg__simple_anndata_with_corpora_and_sparse(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 100)

        self._validate_cxg_and_h5ad_content_match(self.sample_h5ad_filename, self.sample_output_directory, True)

    def test__to_cxg__simple_anndata_no_corpora_and_dense(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 0)

        self._validate_cxg_and_h5ad_content_match(self.sample_h5ad_filename, self.sample_output_directory, False)

    def test__to_cxg__simple_anndata_with_corpora_and_dense(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 0)

        self._validate_cxg_and_h5ad_content_match(self.sample_h5ad_filename, self.sample_output_directory, False)

    def test__to_cxg__simple_anndata_with_corpora_and_dense_using_feature_name_var_index(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename, var_index_column_name="feature_name")
        h5ad_file.to_cxg(self.sample_output_directory, 0)

        self._validate_cxg_and_h5ad_content_match(self.sample_h5ad_filename, self.sample_output_directory, False)
        self._validate_cxg_var_index_column_match(
            self.sample_output_directory,
            "feature_name",
        )

    def test__to_cxg__simple_anndata_with_different_var_index_than_h5ad(self):
        h5ad_file = H5ADDataFile(self.sample_h5ad_filename, var_index_column_name="int_category")
        h5ad_file.to_cxg(self.sample_output_directory, 0)

        self._validate_cxg_var_index_column_match(
            self.sample_output_directory,
            "int_category",
        )

    def test__to_cxg__with_sparse_column_encoding(self):
        anndata = self._create_sample_anndata_dataset()
        anndata.X = np.ones((3, 4))
        sparse_with_column_shift_filename = self._write_anndata_to_file(anndata)

        h5ad_file = H5ADDataFile(sparse_with_column_shift_filename)
        h5ad_file.to_cxg(self.sample_output_directory, 50)

        self._validate_cxg_and_h5ad_content_match(
            sparse_with_column_shift_filename, self.sample_output_directory, False, has_column_encoding=True
        )

        # Clean up
        remove(sparse_with_column_shift_filename)

    def test__slash_in_attribute_name(self):
        # tiledb failure repro code from https://github.com/TileDB-Inc/TileDB-Py/issues/294
        # this fails (throws a TileDBError) for 0.8.0 and below but works for 0.8.1+

        col_name = "fo/o"

        attrs = [tiledb.Attr(name=col_name, dtype=np.int)]
        domain = tiledb.Domain(tiledb.Dim(domain=(0, 99), tile=100, dtype=np.uint32))
        schema = tiledb.ArraySchema(
            domain=domain, sparse=False, attrs=attrs, cell_order="row-major", tile_order="row-major"
        )
        tiledb.DenseArray.create("foo", schema)

        try:
            with tiledb.open("foo", mode="w") as A:
                value = dict()
                value[col_name] = np.zeros((100,), dtype=np.int)
                A[:] = value  # if there's a regression, this statement will throw a TileDBError
                # if we get here we're good
        finally:
            rmtree("foo")

    def _validate_cxg_and_h5ad_content_match(self, h5ad_filename, cxg_directory, is_sparse, has_column_encoding=False):
        anndata_object = anndata.read_h5ad(h5ad_filename)

        # Array locations
        metadata_array_location = f"{cxg_directory}/cxg_group_metadata"
        main_x_array_location = f"{cxg_directory}/X"
        main_xr_array_location = f"{cxg_directory}/Xr"
        main_xc_array_location = f"{cxg_directory}/Xc"
        embedding_array_location = f"{cxg_directory}/emb"
        specific_embedding_array_location = f"{self.sample_output_directory}/emb/awesome_embedding"
        obs_array_location = f"{cxg_directory}/obs"
        var_array_location = f"{cxg_directory}/var"

        # Assert CXG structure
        self.assertEqual(tiledb.object_type(cxg_directory), "group")
        self.assertEqual(tiledb.object_type(obs_array_location), "array")
        self.assertEqual(tiledb.object_type(var_array_location), "array")
        if is_sparse:
            self.assertEqual(tiledb.object_type(main_xr_array_location), "array")
            self.assertEqual(tiledb.object_type(main_xc_array_location), "array")
        else:
            self.assertEqual(tiledb.object_type(main_x_array_location), "array")
        self.assertEqual(tiledb.object_type(embedding_array_location), "group")
        self.assertEqual(tiledb.object_type(specific_embedding_array_location), "array")

        # Validate metadata
        with tiledb.open(metadata_array_location, mode="r") as metadata_array:
            self.assertIn("cxg_version", metadata_array.meta)

        # Validate obs index
        with tiledb.open(obs_array_location, mode="r") as obs_array:
            expected_index_data = anndata_object.obs.index.to_numpy()
            obs_array_schema = json.loads(obs_array.meta["cxg_schema"])
            index_name = obs_array_schema["index"]
            actual_index_data = obs_array.query(attrs=[index_name])[:][index_name]
            self.assertTrue(np.array_equal(expected_index_data, actual_index_data))

            # Validate obs columns
            expected_columns = list(anndata_object.obs.columns.values)
            for column_name in expected_columns:
                if obs_array_schema.get(column_name, {}).get("type", "") == "categorical":
                    categories = obs_array_schema[column_name]["categories"]
                    mapping = dict(zip(range(len(categories)), categories))
                    # Validate values
                    expected_data = anndata_object.obs[column_name].to_numpy()
                    actual_data = obs_array.query(attrs=[column_name])[:][column_name]
                    actual_data = np.array([mapping[i] for i in actual_data], dtype=expected_data.dtype)
                    self.assertTrue(np.array_equal(expected_data, actual_data))

                    # Validate codes
                    expected_data = Categorical(anndata_object.obs[column_name].to_numpy()).codes
                    actual_data = obs_array.query(attrs=[column_name])[:][column_name]
                    self.assertTrue(np.array_equal(expected_data, actual_data))
                else:
                    expected_data = anndata_object.obs[column_name].to_numpy()
                    actual_data = obs_array.query(attrs=[column_name])[:][column_name]
                    self.assertTrue(np.array_equal(expected_data, actual_data))

        # Validate var index
        with tiledb.open(var_array_location, mode="r") as var_array:
            expected_index_data = anndata_object.var.index.to_numpy()
            index_name = json.loads(var_array.meta["cxg_schema"])["index"]
            actual_index_data = var_array.query(attrs=[index_name])[:][index_name]
            self.assertTrue(np.array_equal(expected_index_data, actual_index_data))

            # Validate var columns
            expected_columns = anndata_object.var.columns.values
            for column_name in expected_columns:
                expected_data = anndata_object.var[column_name].to_numpy()
                actual_data = var_array.query(attrs=[column_name])[:][column_name]
                self.assertTrue(np.array_equal(expected_data, actual_data))

        # Validate embedding
        expected_embedding_data = anndata_object.obsm.get("X_awesome_embedding")
        with tiledb.open(specific_embedding_array_location, mode="r") as embedding_array:
            actual_embedding_data = embedding_array[:, 0:2]
            self.assertTrue(np.array_equal(expected_embedding_data, actual_embedding_data))

        # Validate X matrix if not column shifted
        if not has_column_encoding and not is_sparse:
            expected_x_data = anndata_object.X
            with tiledb.open(main_x_array_location, mode="r") as x_array:
                if is_sparse:
                    actual_x_data = np.zeros_like(expected_x_data)
                    data = x_array[:]
                    actual_x_data[data["obs"], data["var"]] = data[""]
                else:
                    actual_x_data = x_array[:, :]
                self.assertTrue(np.array_equal(expected_x_data, actual_x_data))
        elif not has_column_encoding:
            expected_x_data = anndata_object.X
            with tiledb.open(main_xr_array_location, mode="r") as x_array:
                actual_x_data = np.zeros_like(expected_x_data)
                data = x_array[:]
                actual_x_data[data["obs"], data["var"]] = data[""]
                self.assertTrue(np.array_equal(expected_x_data, actual_x_data))

            with tiledb.open(main_xc_array_location, mode="r") as x_array:
                actual_x_data = np.zeros_like(expected_x_data)
                data = x_array[:]
                actual_x_data[data["obs"], data["var"]] = data[""]
                self.assertTrue(np.array_equal(expected_x_data, actual_x_data))

    def _validate_cxg_var_index_column_match(self, cxg_directory, expected_index_name):
        var_array_location = f"{cxg_directory}/var"
        with tiledb.open(var_array_location, mode="r") as var_array:
            actual_index_name = json.loads(var_array.meta["cxg_schema"])["index"]
            self.assertEqual(actual_index_name, expected_index_name)

    def _write_anndata_to_file(self, anndata):
        temporary_filename = fixture_file_path(f"{uuid4()}.h5ad")
        anndata.write(temporary_filename)

        return temporary_filename

    def _create_sample_anndata_dataset(self):
        # Create X
        X = np.random.rand(3, 4)

        # Create obs
        random_string_category = Series(data=["a", "b", "b"], dtype="category")
        random_float_category = Series(data=[3.2, 1.1, 2.2], dtype=np.float32)
        obs_dataframe = DataFrame(
            data={"string_category": random_string_category, "float_category": random_float_category}
        )
        obs = obs_dataframe

        # Create vars
        random_int_category = Series(data=[3, 1, 2, 4], dtype=np.int32)
        random_bool_category = Series(data=[True, True, False, True], dtype=np.bool_)
        feature_name = Series(data=["a", "b", "c", "d"])
        var_dataframe = DataFrame(
            data={
                "int_category": random_int_category.values,
                "bool_category": random_bool_category.values,
                "feature_name": feature_name.values,
            },
            index=feature_name,
        )
        var_dataframe.index.name = "feature_name"
        var = var_dataframe

        # Create embeddings
        random_embedding = np.random.rand(3, 2)
        obsm = {"X_awesome_embedding": random_embedding}

        # Create uns corpora metadata
        uns = {}
        for metadata_field in CorporaConstants.REQUIRED_SIMPLE_METADATA_FIELDS:
            uns[metadata_field] = "random"

        uns["batch_condition"] = np.array(["a", "b"], dtype="object")

        # Need to carefully set the corpora schema versions in order for tests to pass.
        uns["schema_version"] = "4.0.0"

        return anndata.AnnData(X=X, obs=obs, var=var, obsm=obsm, uns=uns)
