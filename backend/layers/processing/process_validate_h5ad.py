from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    CollectionVersionId,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersionId,
)
from backend.layers.common.ingestion_manifest import IngestionManifest
from backend.layers.processing.exceptions import ValidationFailed
from backend.layers.processing.logger import logit
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.thirdparty.s3_provider import S3ProviderInterface
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface


class ProcessValidateH5AD(ProcessingLogic):
    """
    Base class for handling the `Validate` step of the step function.
    This will:
        1. Download the h5ad artifact
        2. upload the original file to S3
        3. Set DatasetStatusKey.H5AD DatasetValidationStatus.VALIDATING
        4. Validate the h5ad
        5. Set DatasetStatusKey.H5AD DatasetValidationStatus.VALID
        6. Set the DatasetStatusKey.RDS DatasetConversionStatus.SKIPPED
    """

    schema_validator: SchemaValidatorProviderInterface

    def __init__(
        self,
        business_logic: BusinessLogicInterface,
        uri_provider: UriProviderInterface,
        s3_provider: S3ProviderInterface,
        schema_validator: SchemaValidatorProviderInterface,
    ) -> None:
        super().__init__()
        self.business_logic = business_logic
        self.uri_provider = uri_provider
        self.s3_provider = s3_provider
        self.schema_validator = schema_validator

    def upload_raw_h5ad(
        self, dataset_version_id: DatasetVersionId, dataset_uri: str, artifact_bucket: str, key_prefix: str
    ) -> str:
        """
        Upload raw h5ad from dataset_uri to artifact bucket

        :param dataset_version_id:
        :param dataset_uri:
        :param artifact_bucket:
        :param key_prefix:
        :return: local_filename: Local filepath to raw h5ad
        """
        self.update_processing_status(dataset_version_id, DatasetStatusKey.PROCESSING, DatasetProcessingStatus.PENDING)

        # Download the original dataset from Dropbox
        local_filename = self.download_from_source_uri(
            source_uri=dataset_uri,
            local_path=CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME,
        )

        # Upload the original dataset to the artifact bucket
        self.update_processing_status(dataset_version_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)
        self.create_artifact(
            local_filename,
            DatasetArtifactType.RAW_H5AD,
            key_prefix,
            dataset_version_id,
            artifact_bucket,
            DatasetStatusKey.H5AD,
        )
        self.update_processing_status(dataset_version_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADED)

        return local_filename

    @logit
    def validate_h5ad_file(self, dataset_version_id: DatasetVersionId, local_filename: str) -> None:
        """
        Validates the specified dataset file and updates the processing status in the database
        :param dataset_version_id: version ID of the dataset to update
        :param local_filename: file name of the dataset to validate and label
        :return: boolean indicating if seurat conversion is possible
        """
        self.update_processing_status(
            dataset_version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALIDATING
        )

        try:
            is_valid, errors, can_convert_to_seurat = self.schema_validator.validate_anndata(local_filename)
        except Exception as e:
            self.logger.exception("validation failed")
            self.update_processing_status(
                dataset_version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.INVALID
            )
            raise ValidationFailed([str(e)]) from None

        if not is_valid:
            self.update_processing_status(
                dataset_version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.INVALID
            )
            raise ValidationFailed(errors)
        else:
            # Skip seurat conversion
            self.update_processing_status(dataset_version_id, DatasetStatusKey.RDS, DatasetConversionStatus.SKIPPED)
            self.update_processing_status(dataset_version_id, DatasetStatusKey.H5AD, DatasetValidationStatus.CONVERTING)

    def process(
        self,
        dataset_version_id: DatasetVersionId,
        manifest: IngestionManifest,
        artifact_bucket: str,
    ):
        """
        1. Download the original dataset from URI
        2. Validate

        :param manifest:
        :param dataset_version_id:
        :param artifact_bucket:
        :return:
        """
        dataset_uri = str(manifest.anndata)
        # validate and upload raw h5ad file to s3
        key_prefix = self.get_key_prefix(dataset_version_id.id)
        local_filename = self.upload_raw_h5ad(dataset_version_id, dataset_uri, artifact_bucket, key_prefix)

        # Validate and label the dataset
        self.validate_h5ad_file(dataset_version_id, local_filename)


class ProcessValidateATACSEQ(ProcessingLogic):
    def __init__(
        self,
        business_logic: BusinessLogicInterface,
        uri_provider: UriProviderInterface,
        s3_provider: S3ProviderInterface,
        schema_validator: SchemaValidatorProviderInterface,
    ) -> None:
        super().__init__()
        self.business_logic = business_logic
        self.uri_provider = uri_provider
        self.s3_provider = s3_provider
        self.schema_validator = schema_validator

    def process(
        self,
        collection_version_id: CollectionVersionId,
        dataset_version_id: DatasetVersionId,
        manifest: IngestionManifest,
        datasets_bucket: str,
    ):
        """

        :param collection_version_id:
        :param dataset_version_id:
        :param manifest:
        :param datasets_bucket:
        :return:
        """
        pass
