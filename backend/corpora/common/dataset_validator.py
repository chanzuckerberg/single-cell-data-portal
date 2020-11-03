import logging
import os
import time

import anndata
import s3fs
from numpy import ndarray
from pandas import DataFrame, Series
import typing
from scipy.sparse import spmatrix

from .utils.corpora_constants import CorporaConstants
from .utils.math_utils import sizeof_formatted


class DatasetValidator:
    """Validates a dataset file that has been uploaded by a submitted to ensure that the correct required metadata
    has been inputted in the expected locations of the file (based on file type) and ensures that no PII exists in
    the dataset file."""

    def __init__(self, s3_uri):
        self.s3_uri = s3_uri
        self.s3_path = s3_uri.replace("s3://", "")
        self.s3_file_system = s3fs.S3FileSystem(
            anon=False, client_kwargs={"endpoint_url": os.getenv("BOTO_ENDPOINT_URL")}
        )

        # Read in ontologies
        for ontology in CorporaConstants.CORPORA_ONTOLOGIES:
            logging.info(f"Reading in {ontology.ontology_name.upper()} ontology.")
            start_time = time.time()

            ontology_file_object = self.s3_file_system.open(ontology.s3_uri, "r")
            ontology_terms_list = ontology_file_object.read().split("\n")
            setattr(self, ontology.ontology_name, ontology_terms_list)
            ontology_file_object.close()

            logging.info(
                f"Completed reading {len(getattr(self, ontology.ontology_name))} values for {ontology.ontology_name} "
                f"ontology in {time.time() - start_time:.3f} seconds."
            )

    def validate_dataset_file(self, loom_x_layer_name=None):
        """
        Reads in file object and triages for specific file type validation.
        """

        file_object = self.s3_file_system.open(self.s3_path, "rb")
        file_object_size = file_object.info().get("Size")
        logging.info(f"Validating file {self.s3_uri} with size {sizeof_formatted(file_object_size)}")

        if self.s3_path.endswith(CorporaConstants.H5AD_FILE_TYPE):
            self.validate_h5ad_dataset(file_object)

        elif self.s3_path.endswith(CorporaConstants.LOOM_FILE_TYPE):
            self.validate_loom_dataset(file_object, loom_x_layer_name)

        else:
            logging.warning(f"Unknown type of dataset with path {self.s3_path}!")

        file_object.close()

    def validate_h5ad_dataset(self, file_object):
        """
        Reads the H5AD file contents into an AnnData object. Each attribute of the AnnData object will then be
        checked to ensure it contains the appropriate metadata.
        """

        start_time = time.time()
        logging.info("Reading H5AD file into anndata object...")
        anndata_object = anndata.read_h5ad(file_object)
        logging.info(f"Finished reading anndata object in {time.time() - start_time:.3f} seconds.")

        self.validate_anndata_object(anndata_object)

    def validate_loom_dataset(self, file_object, loom_x_layer_name=None):
        """
        Reads the Loom file contents into an AnnData object. Each attribute of the AnnData object will then be
        checked to ensure it contains the appropriate metadata.
        """

        start_time = time.time()
        logging.info("Reading Loom file into anndata object...")
        if loom_x_layer_name:
            anndata_object = anndata.read_loom(file_object, X_name=loom_x_layer_name)
        else:
            anndata_object = anndata.read_loom(file_object)
        logging.info(f"Finished reading anndata object in {time.time() - start_time:.3f} seconds.")

        self.validate_anndata_object(anndata_object)

    def validate_anndata_object(self, anndata_object: anndata.AnnData):
        start_time = time.time()
        logging.info("Beginning validation of anndata object...")
        self.verify_layers(anndata_object)
        self.verify_obs(anndata_object)
        self.verify_vars(anndata_object)
        self.verify_uns(anndata_object)
        logging.info(f"Finished completing validation in {time.time() - start_time:.3f} seconds.")

    def verify_layers(self, data_object: anndata.AnnData):
        """
        Verifies that the dataset contains at least the raw data and if other layers are provided, that they each
        contain an appropriate description.
        """

        # Check to make sure X data exists
        has_data = True
        if isinstance(data_object.X, DataFrame):
            has_data = data_object.X.data.any()
        elif isinstance(data_object.X, ndarray):
            has_data = data_object.X.any()
        elif isinstance(data_object.X, spmatrix):
            has_data = (data_object.X.count_nonzero() == data_object.X.nnz) or data_object.X.nnz == 0
        else:
            logging.warning(
                f"Could not check X data layer to ensure that it exists. The type is " f"{type(data_object.X)}!"
            )

        if not has_data:
            logging.warning("No data in the X layer can be found in the dataset or all observations are zeros!")

        # Ensure that the layer_descriptions metadata key exists in the `uns` field of the anndata object.
        if (CorporaConstants.LAYER_DESCRIPTIONS not in data_object.uns_keys()) or (
            not data_object.uns.get(CorporaConstants.LAYER_DESCRIPTIONS)
        ):
            logging.warning("Required layers descriptions are missing from uns field to describe data layers!")
        else:
            # Check to ensure that there are descriptions for each layer
            for layer_name in data_object.layers.keys():
                if layer_name not in data_object.uns.get(CorporaConstants.LAYER_DESCRIPTIONS).keys():
                    logging.warning(f"Missing layer description for layer {layer_name}!")

            # Check to make sure that X has a layer description and if the anndata populate the `raw` field,
            # that a raw data layer description also exists.
            if (
                CorporaConstants.X_DATA_LAYER_NAME
                not in data_object.uns.get(CorporaConstants.LAYER_DESCRIPTIONS).keys()
            ):
                logging.warning(f"Missing layer description for layer {CorporaConstants.X_DATA_LAYER_NAME}!")
            if data_object.raw:
                if (
                    CorporaConstants.RAW_DATA_LAYER_NAME
                    not in data_object.uns.get(CorporaConstants.LAYER_DESCRIPTIONS).keys()
                ):
                    logging.warning(f"Missing layer description for layer {CorporaConstants.RAW_DATA_LAYER_NAME}!")

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
            if metadata_field.field_name not in observation_keys:
                self.log_error_message(metadata_field.field_name, "obs", type(data_object).__name__)
            else:
                self.verify_metadata_type(metadata_field, data_object.obs.get(metadata_field.field_name))

    def verify_vars(self, data_object: anndata.AnnData):
        """
        Validates the variable attribute of the AnnData object to ensure that all variable IDs are unique.
        """

        if data_object.var.index.duplicated().any():
            logging.warning("Each variable is not unique!")

    def verify_uns(self, data_object: anndata.AnnData):
        """
        Validate the unstructured attribute of the AnnData object to ensure that it contains the appropriate
        dataset-level and collection-level metadata and outputs which metadata fields are missing. Note that no
        exception is thrown when metadata is found to be missing and rather an informative message is outputted instead.
        """

        unstructured_metadata_keys = data_object.uns_keys()

        for metadata_field in (
            CorporaConstants.REQUIRED_DATASET_METADATA_FIELDS
            + CorporaConstants.REQUIRED_DATASET_PRESENTATION_METADATA_FIELDS
        ):
            if metadata_field.field_name not in unstructured_metadata_keys:
                self.log_error_message(metadata_field.field_name, "uns", type(data_object).__name__)
            else:
                self.verify_metadata_type(metadata_field, data_object.uns.get(metadata_field.field_name))

    def verify_metadata_type(
        self,
        metadata_property: CorporaConstants.TypedMetadata,
        metadata_values_in_dataset: typing.Union[Series, str, list, dict],
    ):
        """
        Validates the type of each value in a property of an AnnData object.

        Each property value passed in the pandas Series object `metadata_values_in_dataset` is expected to be of the
        type `metadata_property` where `metadata_property` can either be a type or can be a namedtuple that represents
        an ontology. When the type is an ontology, the validation instead ensures that each value of the AnnData
        property belongs to the expected ontology. In some cases, there is an acceptable alternate value not part of
        the ontology for which the validation also checks.
        """

        # Canonicalize the metadata values type into a list
        if isinstance(metadata_values_in_dataset, Series):
            metadata_values = metadata_values_in_dataset.values
        else:
            metadata_values = [metadata_values_in_dataset]

        # Handle the case where the required type is simply a type check.
        if isinstance(metadata_property.required_type, type):
            for data_value in metadata_values:
                if not isinstance(data_value, metadata_property.required_type):
                    logging.warning(
                        f"Value {data_value} of type {type(data_value)} is not of expected type "
                        f"{metadata_property.required_type} for metadata field {metadata_property.field_name}."
                    )

        # Handle the case where the required type is an Ontology.
        elif isinstance(metadata_property.required_type, CorporaConstants.Ontology):
            valid_ontology_names = getattr(self, metadata_property.required_type.ontology_name)
            unrecognized_data_values = set()
            for data_value in metadata_values:
                if not (
                    data_value in valid_ontology_names
                    or (
                        metadata_property.valid_alternative is not None
                        and data_value == metadata_property.valid_alternative
                    )
                ):
                    unrecognized_data_values.add(data_value)

            if unrecognized_data_values:
                logging.warning(
                    f"Values {unrecognized_data_values} were not recognized as a valid value in the "
                    f"{metadata_property.required_type.ontology_name} ontology."
                )

        # Handle the case where the required type is an enum which is represented by a list of accepted values.
        elif isinstance(metadata_property.required_type, list):
            unrecognized_data_values = set()
            for data_value in metadata_values:
                if data_value not in metadata_property.required_type:
                    unrecognized_data_values.add(data_value)

            if unrecognized_data_values:
                logging.warning(
                    f"Values {unrecognized_data_values} are not part of the accepted enum values "
                    f"{metadata_property.required_type} for metadata field {metadata_property.field_name}."
                )

        else:
            logging.warning(
                f"Unable to parse metadata property: {metadata_property.field_name} with type "
                f"{type(metadata_property.required_type)}"
            )

    def log_error_message(self, metadata_field_name, expected_location, dataset_type):
        """
        Pretty-printer of missing metadata fields errors.
        """

        is_ontology = " ontology " if "ONTOLOGY" in metadata_field_name else " "
        logging.warning(
            f"ERROR: Missing{is_ontology}metadata field {metadata_field_name} from {expected_location} in "
            f"{dataset_type} file!"
        )
