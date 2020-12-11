import random
import string
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

        self._create_fake_ontologies_in_s3_bucket()

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
    def test__validate_h5ad_dataset__empty_x_data_outputs_error(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object(0, 0)

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "empty.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("No data in the X layer can be found", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__all_zeros_x_data_outputs_error(self, mock_read_anndata):
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
    def test__validate_h5ad_dataset__validate_different_types_of_x_data_passes(
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
    def test__validate_h5ad_dataset__unsupported_x_data_type_outputs_error(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        test_anndata.X = test_anndata.X.tostring()

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "empty.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("Could not check X data layer", logger.output[0])

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

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__incorrect_metadata_value_type(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        # Expected type of tissue is string so inputting integers should cause an error to be outputted.
        test_anndata.obs["tissue"] = [random.randint(1, 10) for _ in range(len(test_anndata.obs.get("tissue")))]

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "metadata_values.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("is not of expected type", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__unaccepted_ontology_value(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        # Expected type of tissue ontology is to be part of UBERON ontology which is generated in this test suite to be
        # 8 random characters long. Therefore inputting "unknown" which is 7 letters, should cause an error.
        test_anndata.obs["tissue_ontology"] = ["unknown" for _ in range(len(test_anndata.obs.get("tissue_ontology")))]

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "metadata_values.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn(
                f"not recognized as a valid value in the {CorporaConstants.UBERON_ONTOLOGY.ontology_name} ontology",
                logger.output[0],
            )

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__unaccepted_enum_value(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        # Expected enums for "sex" do not include asdfglkjh, so setting the attribute to that should throw an error.
        test_anndata.obs["sex"] = ["asdfglkjh" for _ in range(len(test_anndata.obs.get("sex")))]

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "metadata_values.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("not part of the accepted enum values", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__unaccepted_uns_metadata_type(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        # Expected type of the "title" is string so inputting an integer should cause an error.
        test_anndata.uns["title"] = 12345

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "metadata_values.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn("is not of expected type", logger.output[0])

    @patch("anndata.read_h5ad")
    def test__validate_h5ad_dataset__unaccepted_uns_ontology(self, mock_read_anndata):
        test_anndata = self._create_fully_populated_anndata_object()
        # Expected type of organism ontology is to be part of NCBI Taxon ontology which is generated in this test suite
        # to be 8 random characters long. Therefore inputting "unknown" which is 7 letters, should cause an error.
        test_anndata.uns["organism_ontology"] = "unknown"

        mock_read_anndata.return_value = test_anndata

        dataset_filename = "metadata_values.h5ad"
        s3_uri = self._create_dataset_object_in_s3_bucket(dataset_filename)

        # Run validator
        validator = DatasetValidator(s3_uri)

        # Validate result
        with self.assertLogs(level="WARN") as logger:
            validator.validate_dataset_file()
            self.assertIn(
                f"not recognized as a valid value in the {CorporaConstants.NCBI_TAXON_ONTOLOGY.ontology_name} ontology",
                logger.output[0],
            )

    def _create_fully_populated_anndata_object(self, obs_count=3, var_count=4):
        """
        Create an anndata with `obs_count` observations (usually cells) and `var_count` variables (usually genes).
        Populate every required piece of metadata by Corpora into the anndata as well in the `obs` and `uns` views
        filled with bogus values.
        """

        # Format metadata for obs field of anndata object and fill in with bogus values
        obs_data = {}
        for metadata_field in (
            CorporaConstants.REQUIRED_OBSERVATION_METADATA_FIELDS
            + CorporaConstants.REQUIRED_OBSERVATION_ONTOLOGY_METADATA_FIELDS
        ):
            if metadata_field.required_type is str:
                obs_data[metadata_field.field_name] = [
                    "".join(random.sample(string.ascii_lowercase, 8)) for _ in range(obs_count)
                ]
            elif type(metadata_field.required_type) is list:
                obs_data[metadata_field.field_name] = random.sample(metadata_field.required_type, obs_count)
            elif isinstance(metadata_field.required_type, CorporaConstants.Ontology):
                obs_data[metadata_field.field_name] = random.sample(
                    getattr(self, metadata_field.required_type.ontology_name), obs_count
                )

        # Format unstructured metadata for uns field of anndata object and fill in with bogus values
        uns_data = {}
        for metadata_field in (
            CorporaConstants.REQUIRED_DATASET_PRESENTATION_METADATA_FIELDS
            + CorporaConstants.REQUIRED_DATASET_PRESENTATION_HINTS_METADATA_FIELDS
            + CorporaConstants.REQUIRED_DATASET_METADATA_FIELDS
            + CorporaConstants.OPTIONAL_COLLECTION_LEVEL_METADATA_FIELDS
        ):
            if metadata_field.required_type is str:
                uns_data[metadata_field.field_name] = "".join(random.sample(string.ascii_lowercase, 8))
            elif metadata_field.required_type is list:
                uns_data[metadata_field.field_name] = []
            elif metadata_field.required_type is dict:
                uns_data[metadata_field.field_name] = {}
            elif isinstance(metadata_field.required_type, CorporaConstants.Ontology):
                uns_data[metadata_field.field_name] = random.choice(
                    getattr(self, metadata_field.required_type.ontology_name)
                )

        # Add layers descriptions
        uns_data[CorporaConstants.LAYER_DESCRIPTIONS].update({CorporaConstants.X_DATA_LAYER_NAME: "random layer"})

        # Generate random data
        data = numpy.random.rand(obs_count, var_count)

        return anndata.AnnData(X=data, obs=DataFrame(data=obs_data), uns=uns_data)

    def _create_fake_ontologies_in_s3_bucket(self):
        """
        Create fake ontologies in the bogus S3 resource.
        """

        for ontology in CorporaConstants.CORPORA_ONTOLOGIES:
            number_of_ontology_terms = random.randint(5, 10)
            ontology_term_list = [
                "".join(random.sample(string.ascii_lowercase, 8)) for _ in range(number_of_ontology_terms)
            ]
            setattr(self, ontology.ontology_name, ontology_term_list)
            self._create_dataset_object_in_s3_bucket(
                ontology.s3_uri.split("/")[1], ontology.s3_uri.split("/")[0], "\n".join(ontology_term_list)
            )

    def _create_dataset_object_in_s3_bucket(self, filename, bucket_name=None, content="file_contents"):
        dataset_s3_object = self.create_s3_object(filename, bucket_name=bucket_name, content=content)
        return "s3://" + dataset_s3_object.bucket_name + "/" + dataset_s3_object.key
