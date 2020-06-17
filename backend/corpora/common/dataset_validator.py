import logging
import time

import anndata
import s3fs

from .utils.corpora_constants import CorporaConstants
from .utils.math_utils import sizeof_formatted


class DatasetValidator:
    """ Validates a dataset file that has been uploaded by a submitted to ensure that the correct required metadata
    has been inputted in the expected locations of the file (based on file type) and ensures that no PII exists in
    the dataset file."""

    def __init__(self, s3_uri):
        self.s3_uri = s3_uri
        self.s3_path = s3_uri.replace("s3://", "")
        self.s3_file_system = s3fs.S3FileSystem()

    def validate_dataset_file(self):
        """
        Reads in file object and triages for specific file type validation.
        """

        file_object = self.s3_file_system.open(self.s3_path, "rb")
        file_object_size = file_object.info().get("Size")
        logging.info(f"Validating file {self.s3_uri} with size {sizeof_formatted(file_object_size)}")

        if self.s3_path.endswith(CorporaConstants.H5AD_FILE_TYPE):
            self.validate_h5ad_dataset(file_object)

        elif self.s3_path.endswith(CorporaConstants.LOOM_FILE_TYPE):
            self.validate_loom_dataset(file_object)

        else:
            logging.warning(f"Unknown type of dataset with path {self.s3_path}!")

        file_object.close()

    def validate_h5ad_dataset(self, file_object):
        """
        Reads the file contents into an AnnData object. Each attribute of the AnnData object will then be checked to
        ensure it contains the appropriate metadata.
        """

        start_time = time.time()
        logging.info("Reading anndata object...")
        anndata_object = anndata.read_h5ad(file_object)
        logging.info(f"Finished reading anndata object in {time.time() - start_time:.3f} seconds.")

        start_time = time.time()
        logging.info("Beginning validation of file...")
        self.verify_obs(anndata_object)
        self.verify_vars(anndata_object)
        self.verify_uns(anndata_object)
        logging.info(f"Finished completing validation on the file in {time.time() - start_time:.3f} seconds.")

    def verify_obs(self, data_object: anndata.AnnData):
        """
        Validates the observation attribute of an AnnData object. Checks to ensure that all observation IDs are
        unique and that the observation metadata fields as described by the Corpora Schema exist. If the validation
        fails in any way, the errors are outputted rather than the validation aborted.
        """

        observation_keys = data_object.obs_keys()

        # Check to ensure that all IDs are unique
        if data_object.obs.index.duplicated().any():
            logging.warning("Each observation is not unique!")

        for metadata_field in (
            CorporaConstants.REQUIRED_OBSERVATION_METADATA_FIELDS
            + CorporaConstants.REQUIRED_OBSERVATION_ONTOLOGY_METADATA_FIELDS
        ):
            if metadata_field not in observation_keys:
                self.log_error_message(metadata_field, "obs", type(data_object).__name__)

    def verify_vars(self, data_object: anndata.AnnData):
        """
        Validates the variable attribute of the AnnData object to ensure that all variable IDs are unique.
        """

        if data_object.var.index.duplicated().any():
            logging.warning("Each variable is not unique!")

    def verify_uns(self, data_object: anndata.AnnData):
        """
        Validate the unstructured attribute of the AnnData object to ensure that it contains the appropriate
        dataset-level and project-level metadata and outputs which metadata fields are missing. Note that no
        exception is thrown when metadata is found to be missing and rather an informative message is outputted instead.
        """

        unstructured_metadata_keys = data_object.uns_keys()

        for metadata_field in (
            CorporaConstants.REQUIRED_DATASET_METADATA_FIELDS
            + CorporaConstants.REQUIRED_DATASET_PRESENTATION_METADATA_FIELDS
        ):
            if metadata_field not in unstructured_metadata_keys:
                self.log_error_message(metadata_field, "uns", type(data_object).__name__)

    def log_error_message(self, metadata_field_name, expected_location, dataset_type):
        """
        Pretty-printer of missing metadata fields errors.
        """

        is_ontology = " ontology " if "ONTOLOGY" in metadata_field_name else " "
        logging.warning(
            f"ERROR: Missing{is_ontology}metadata field {metadata_field_name} from {expected_location} in "
            f"{dataset_type} file!"
        )

    def validate_loom_dataset(self, file_object):
        # TODO: Implement this as part of ticket #375.
        pass
