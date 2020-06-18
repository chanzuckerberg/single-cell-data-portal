import random
from unittest.mock import patch

import anndata
import numpy
from pandas import DataFrame
from scipy.sparse.csr import csr_matrix

from backend.corpora.common.dataset_validator import DatasetValidator
from backend.corpora.common.utils.corpora_constants import CorporaConstants
from .. import CorporaTestCaseUsingMockAWS


class TestDatasetValidator(CorporaTestCaseUsingMockAWS):
    def setUp(self):
        super().setUp()

    def test__validate_dataset_file__unknown_file(self):
        # Setup mocks
        unaccepted_dataset_filename = "dont_process_me.pdf"
        unaccepted_s3_uri = self._create_dataset_object_in_s3_bucket(unaccepted_dataset_filename)

        # Run validator
        validator = DatasetValidator(unaccepted_s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("Unknown type of dataset", logger.output[0])

    @patch("logging.warning")
    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__contains_all_metadata(self, mock_read_anndata, mock_log_warning):
        test_anndata = self._create_fully_populated_anndata_object()

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "fully_populated_h5ad.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)
        validator.validate_dataset_file()

        # Validate result
        assert not mock_log_warning.called

    @patch("logging.warning")
    @patch("anndata.read_loom")
    def test__validate_loom_dataset__contains_all_metadata(self, mock_read_anndata, mock_log_warning):
        test_anndata = self._create_fully_populated_anndata_object()

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "fully_populated_loom.loom"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)
        validator.validate_dataset_file()

        # Validate result
        assert not mock_log_warning.called

    @patch("logging.warning")
    @patch("anndata.read_loom")
    def test__validate_loom_dataset__specific_x_layer(self, mock_read_anndata, mock_log_warning):
        test_anndata = self._create_fully_populated_anndata_object()

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "fully_populated_loom.loom"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)
        validator.validate_dataset_file(loom_x_layer_name="main_x")

        # Validate result
        assert not mock_log_warning.called

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__test_case_sensitivity_outputs_error(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()

        # Capitalize one piece of metadata which will output an error
        missing_metadata = "tissue"
        del test_anndata.obs[missing_metadata]
        test_anndata.obs["TISSUE"] = None

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "missing_one_metadata_h5ad.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn(f"Missing metadata field {missing_metadata} from obs", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__missing_obs_metadata_outputs_error(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()

        # Delete one piece of metadata which will be the missing one
        missing_metadata = "tissue"
        del test_anndata.obs[missing_metadata]

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "missing_one_metadata_h5ad.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn(f"Missing metadata field {missing_metadata} from obs", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__missing_uns_metadata_outputs_error(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()

        # Delete one piece of metadata which will be the missing one
        missing_metadata = "organism"
        del test_anndata.uns[missing_metadata]

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "missing_one_metadata_h5ad.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn(f"Missing metadata field {missing_metadata} from uns", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__non_unique_var_outputs_error(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()

        # Create data frame for vars that contains non unique indices
        data_frame = DataFrame(index=["0", "1", "2", "2"])

        test_anndata.var = data_frame

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "non_unique_var.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("Each variable is not unique", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__non_unique_obs_outputs_error(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()

        # Create data frame for vars that contains non unique indices
        data_frame = DataFrame(index=["0", "1", "1"])

        test_anndata.obs = data_frame

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "non_unique_obs.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("Each observation is not unique", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__empty_raw_data_outputs_error(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object(0, 0)

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "empty.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("No raw data can be found", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__all_zeros_raw_data_outputs_error(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        for i in range(test_anndata.X.shape[0]):
            for j in range(test_anndata.X.shape[1]):
                test_anndata.X[i][j] = 0

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "empty.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("all observations are zeros", logger.output[0])

    @patch("logging.warning")
    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__validate_different_types_of_raw_data_passes(
        self, mock_read_anndata, mock_log_warning
    ):
        dataset_filename = "data.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Test with ndarray type
        test_anndata = self._create_fully_populated_anndata_object()
        ndarray_formatted_raw_data = test_anndata.X
        mock_read_anndata.return_value = test_anndata

        validator = DatasetValidator(s3_uri)
        validator.validate_dataset_file()
        assert not mock_log_warning.called

        # Test with DataFrame
        test_anndata.X = DataFrame(data=ndarray_formatted_raw_data)
        mock_read_anndata.return_value = test_anndata
        validator.validate_dataset_file()
        assert not mock_log_warning.called

        # Test with csr_matrix
        test_anndata.X = csr_matrix(ndarray_formatted_raw_data)
        mock_read_anndata.return_value = test_anndata
        validator.validate_dataset_file()
        assert not mock_log_warning.called

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__missing_layer_description(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        del test_anndata.uns[CorporaConstants.LAYER_DESCRIPTIONS]

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "layers.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("Required layers descriptions are missing", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__missing_one_specific_layer_description(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        fake_layer_one = numpy.random.permutation(test_anndata.X)
        fake_layer_two = numpy.random.permutation(fake_layer_one)
        test_anndata.layers["my_layer_one"] = fake_layer_one
        test_anndata.layers["my_layer_two"] = fake_layer_two
        test_anndata.uns[CorporaConstants.LAYER_DESCRIPTIONS].update({"my_layer_one": "my one lovely layer"})

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "layers.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("Missing layer description for layer my_layer_two", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__missing_x_layer_description(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        fake_layer_one = numpy.random.permutation(test_anndata.X)
        test_anndata.layers["my_layer_one"] = fake_layer_one
        test_anndata.uns[CorporaConstants.LAYER_DESCRIPTIONS].update({"my_layer_one": "my one lovely layer"})
        del test_anndata.uns[CorporaConstants.LAYER_DESCRIPTIONS][CorporaConstants.X_DATA_LAYER_NAME]

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "layers.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("Missing layer description for layer X", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__missing_raw_data_layer_description(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        test_anndata.raw = test_anndata

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "layers.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn(
                f"Missing layer description for layer {CorporaConstants.RAW_DATA_LAYER_NAME}", logger.output[0]
            )

    def _create_fully_populated_anndata_object(self, obs_count=3, var_count=4):
        """
        Create an anndata with `obs_count` observations (usually cells) and `var_count` variables (usually genes).
        Populate every required piece of metadata by Corpora into the anndata as well in the `obs` and `uns` views
        filled with bogus values.
        """

        # Format metadata for obs field of anndata object
        obs_data = {}
        for metadata_field in (
            CorporaConstants.REQUIRED_OBSERVATION_METADATA_FIELDS
            + CorporaConstants.REQUIRED_OBSERVATION_ONTOLOGY_METADATA_FIELDS
        ):
            obs_data[metadata_field] = random.sample(range(10, 30), obs_count)

        # Format unstructured metadata for uns field of anndata object
        uns_data = {}
        for metadata_field in (
            CorporaConstants.REQUIRED_DATASET_PRESENTATION_METADATA_FIELDS
            + CorporaConstants.REQUIRED_DATASET_PRESENTATION_HINTS_METADATA_FIELDS
            + CorporaConstants.REQUIRED_DATASET_METADATA_FIELDS
            + CorporaConstants.OPTIONAL_PROJECT_LEVEL_METADATA_FIELDS
        ):
            uns_data[metadata_field] = {}

        # Add layers descriptions
        uns_data[CorporaConstants.LAYER_DESCRIPTIONS][CorporaConstants.X_DATA_LAYER_NAME] = "random layer"

        # Generate random data
        data = numpy.random.rand(obs_count, var_count)

        return anndata.AnnData(X=data, obs=DataFrame(data=obs_data), uns=uns_data)

    def _create_dataset_object_in_s3_bucket(self, filename):
        dataset_s3_object = self.create_s3_object(filename)
        return "s3://" + dataset_s3_object.bucket_name + "/" + dataset_s3_object.key
