import hashlib

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    ARTIFACT_TO_EXTENSION,
    CollectionVersionId,
    DatasetArtifactId,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetStatusKey,
    DatasetValidationStatus,
    DatasetVersionId,
)
from backend.layers.common.ingestion_manifest import IngestionManifest
from backend.layers.processing.exceptions import ConversionFailed, ValidationAtacFailed
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.thirdparty.s3_provider_interface import S3ProviderInterface
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface


class ProcessValidateATAC(ProcessingLogic):
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

    def create_atac_artifact(
        self,
        file_name: str,
        artifact_type: DatasetArtifactType,
        dataset_version_id: DatasetVersionId,
        processing_status_key: DatasetStatusKey,
        datasets_bucket: str,
        fragment_artifact_id: DatasetArtifactId = None,
    ) -> DatasetArtifactId:
        """
        Uploads the file to S3 and updates the database with the artifact using the artifact_id as the prefix for the
        key.
        :param file_name: the local file to upload
        :param artifact_type: the type of artifact to upload
        :param dataset_version_id: the dataset version id
        :param processing_status_key: the key to update the processing status
        :param datasets_bucket: the bucket to upload the dataset to
        :param fragment_artifact_id: the artifact id of the fragment file to be use in the fragment index file.
        :return:
        """
        self.update_processing_status(dataset_version_id, processing_status_key, DatasetConversionStatus.UPLOADING)
        try:
            artifact_id = self.business_logic.add_dataset_artifact(dataset_version_id, artifact_type, "dummy")
            if fragment_artifact_id:
                key_prefix = self.get_key_prefix(fragment_artifact_id.id)
            else:
                key_prefix = self.get_key_prefix(artifact_id.id)
            key = ".".join((key_prefix, ARTIFACT_TO_EXTENSION[artifact_type]))
            self.s3_provider.upload_file(
                file_name, datasets_bucket, key, extra_args={"ACL": "bucket-owner-full-control"}
            )
            datasets_s3_uri = self.make_s3_uri(datasets_bucket, key_prefix, key)
            self.logger.info(f"Uploaded {dataset_version_id}.{artifact_type} to {datasets_s3_uri}")
            self.business_logic.update_dataset_artifact(artifact_id, datasets_s3_uri)
            self.logger.info(f"Updated database with {artifact_type}.")
            self.update_processing_status(dataset_version_id, processing_status_key, DatasetConversionStatus.UPLOADED)
            return artifact_id
        except Exception as e:
            self.logger.error(e)
            raise ConversionFailed(processing_status_key) from None

    def skip_atac_validation(
        self, local_anndata_filename: str, manifest: IngestionManifest, dataset_version_id
    ) -> bool:
        """
        Check if atac validation should be skipped
        :param local_anndata_filename: the local anndata file
        :param manifest: the manifest
        :param dataset_version_id: the dataset version id
        :return: True if the validation should be skipped, False otherwise
        """
        # check if the anndata should have a fragment file
        try:
            result = self.schema_validator.check_anndata_requires_fragment(local_anndata_filename)
        except ValueError as e:  # fragment file forbidden
            self.logger.warning(f"Anndata does not support atac fragment files for the follow reason: {e}")
            self.logger.warning("Skipping fragment validation")
            self.update_processing_status(
                dataset_version_id,
                DatasetStatusKey.ATAC_FRAGMENT,
                DatasetConversionStatus.SKIPPED,
                validation_atac_errors=[str(e)],
            )
            return True

        if manifest.atac_fragment is None:
            if result:  # fragment file required
                self.update_processing_status(
                    dataset_version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.INVALID
                )
                raise ValidationAtacFailed(["Anndata requires fragment file"])
            else:  # fragment file optional
                self.logger.info("Fragment is optional and not present. Skipping fragment validation.")
                self.update_processing_status(
                    dataset_version_id,
                    DatasetStatusKey.ATAC_FRAGMENT,
                    DatasetConversionStatus.SKIPPED,
                    validation_atac_errors=["Fragment is optional and not present."],
                )
                return True
        return False

    def hash_file(self, file_name: str) -> str:
        """
        Hash the file
        :param file_name: the file to hash
        :return: the hash
        """
        with open(file_name, "rb") as f:
            return hashlib.sha256(f.read()).hexdigest()

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
        # Download the original dataset files from URI
        local_anndata_filename = self.download_from_source_uri(
            source_uri=str(manifest.anndata),
            local_path=CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME,
        )

        if self.skip_atac_validation(local_anndata_filename, manifest, dataset_version_id):
            return

        # Download the original fragment file from URI
        local_fragment_filename = self.download_from_source_uri(
            source_uri=str(manifest.atac_fragment), local_path=CorporaConstants.ORIGINAL_ATAC_FRAGMENT_FILENAME
        )

        # Hash the original fragment file. This is used to check if the new fragment file is the same as the old
        # fragment file to avoid uploading the same file multiple times
        if not self.business_logic.is_public_uri(manifest.atac_fragment):
            # new fragment file is from a private source, so we must upload it
            original_fragment_hash = None
        else:
            original_fragment_hash = self.hash_file(local_fragment_filename)

        # Validate the fragment with anndata file
        try:
            errors = self.schema_validator.validate_atac(local_fragment_filename, local_anndata_filename)
        except Exception as e:
            # for unexpected errors, log the exception and raise a ValidationAtacFailed exception
            self.logger.exception("validation failed")
            self.update_processing_status(
                dataset_version_id, DatasetStatusKey.ATAC_FRAGMENT, DatasetConversionStatus.FAILED
            )
            raise ValidationAtacFailed([str(e)]) from None

        if errors:
            # if the validation fails, update the processing status and raise a ValidationAtacFailed exception
            self.update_processing_status(
                dataset_version_id, DatasetStatusKey.ATAC_FRAGMENT, DatasetConversionStatus.FAILED
            )
            raise ValidationAtacFailed(errors)

        # check to see if the new fragments is the same as the old fragment
        # if it is the same, skip the upload and use link the old fragment to the new dataset
        # This case is for when the dataset is reprocessed, or updated
        if original_fragment_hash is not None and original_fragment_hash == self.hash_file(local_fragment_filename):
            # get the artifact id of the old fragment, and add it to the new dataset
            artifact_name = str(manifest.atac_fragment).split("/")[-1]
            artifact = self.persistence.get_artifact_by_uri_suffix(artifact_name)
            self.business_logic.add_artifact_to_dataset_version(dataset_version_id, artifact.id)
            # get the artifact id of the old fragment index, and add it to the new dataset
            artifact = self.persistence.get_artifact_by_uri_suffix(artifact_name + ".tbi")
            self.business_logic.add_artifact_to_dataset_version(dataset_version_id, artifact.id)
        else:
            fragment_artifact_id = self.create_atac_artifact(
                local_fragment_filename,
                DatasetArtifactType.ATAC_FRAGMENT,
                dataset_version_id,
                DatasetStatusKey.ATAC_FRAGMENT,
                datasets_bucket,
            )
            self.create_atac_artifact(
                local_fragment_filename + ".tbi",
                DatasetArtifactType.ATAC_INDEX,
                dataset_version_id,
                DatasetStatusKey.ATAC_FRAGMENT,
                datasets_bucket,
                fragment_artifact_id,
            )
        self.logger.info("Processing completed successfully")
        return
