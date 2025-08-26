import hashlib
import os

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
        datasets_bucket: str,
        fragment_artifact_id: DatasetArtifactId = None,
    ) -> DatasetArtifactId:
        """
        Uploads the file to S3 and updates the database with the artifact using the artifact_id as the prefix for the
        key.
        :param file_name: the local file to upload
        :param artifact_type: the type of artifact to upload
        :param dataset_version_id: the dataset version id
        :param datasets_bucket: the bucket to upload the dataset to
        :param fragment_artifact_id: the artifact id of the fragment file, to be used in the fragment index filepath for storage
        :return:
        """
        self.update_processing_status(dataset_version_id, DatasetStatusKey.ATAC, DatasetConversionStatus.UPLOADING)
        try:
            artifact_id = DatasetArtifactId()
            if fragment_artifact_id:
                key_prefix = self.get_key_prefix(fragment_artifact_id.id)
            else:
                key_prefix = self.get_key_prefix(artifact_id.id)
            key = f"{key_prefix}-fragment.{ARTIFACT_TO_EXTENSION[artifact_type]}"
            datasets_s3_uri = self.upload_artifact(file_name, key, datasets_bucket)
            self.logger.info(f"Uploaded [{dataset_version_id}/{artifact_type}] to {datasets_s3_uri}")
            self.business_logic.add_dataset_artifact(dataset_version_id, artifact_type, datasets_s3_uri, artifact_id)
            self.logger.info(f"Updated database with {artifact_type}.")
            return artifact_id
        except Exception as e:
            self.logger.error(e)
            raise ConversionFailed(
                DatasetStatusKey.ATAC,
            ) from None

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
            self.logger.warning(f"Anndata does not support atac fragment files for the following reason: {e}")
            if manifest.atac_fragment:
                self.update_processing_status(
                    dataset_version_id,
                    DatasetStatusKey.VALIDATION,
                    DatasetValidationStatus.INVALID,
                )
                raise ValidationAtacFailed(errors=[str(e), "Fragment file not allowed for non atac anndata."]) from None
            self.logger.warning("Fragment validation not applicable for dataset assay type.")
            self.update_processing_status(
                dataset_version_id,
                DatasetStatusKey.ATAC,
                DatasetConversionStatus.NA,
            )
            return True

        if manifest.atac_fragment is None:
            if result:  # fragment file required
                self.update_processing_status(
                    dataset_version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.INVALID
                )
                raise ValidationAtacFailed(errors=["Anndata requires fragment file"])
            else:  # fragment file optional
                self.logger.info("Fragment is optional and not present. Skipping fragment validation.")
                self.update_processing_status(
                    dataset_version_id,
                    DatasetStatusKey.ATAC,
                    DatasetConversionStatus.SKIPPED,
                )
                return True
        return False

    def hash_file(self, file_name: str) -> str:
        """
        Hash the file
        :param file_name: the file to hash
        :return: the hash
        """
        hashobj = hashlib.md5()
        buffer = bytearray(2**18)
        view = memoryview(buffer)
        with open(file_name, "rb") as f:
            while chunk := f.readinto(buffer):
                hashobj.update(view[:chunk])
        return hashobj.hexdigest()

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

        try:
            # Download the original fragment file from URI
            local_fragment_filename = self.download_from_source_uri(
                source_uri=str(manifest.atac_fragment), local_path=CorporaConstants.ORIGINAL_ATAC_FRAGMENT_FILENAME
            )
        except Exception as e:
            self.logger.exception(f"Failed to download fragment file from {manifest.atac_fragment}")
            self.update_processing_status(dataset_version_id, DatasetStatusKey.ATAC, DatasetConversionStatus.FAILED)
            raise ValidationAtacFailed(errors=[str(e)]) from None

        # Deduplicate the fragment file
        if manifest.flags and manifest.flags.deduplicate_fragment:
            try:
                self.schema_validator.deduplicate_fragment(local_fragment_filename)
            except Exception as e:
                self.logger.exception(f"Failed to deduplicate fragment file {local_fragment_filename}")
                self.update_processing_status(dataset_version_id, DatasetStatusKey.ATAC, DatasetConversionStatus.FAILED)
                raise ValidationAtacFailed(errors=[str(e)]) from None

        # Validate the fragment with anndata file
        try:
            errors, fragment_index_file, fragment_file = self.schema_validator.validate_atac(
                local_fragment_filename, local_anndata_filename, CorporaConstants.NEW_ATAC_FRAGMENT_FILENAME
            )
        except Exception as e:
            # for unexpected errors, log the exception and raise a ValidationAtacFailed exception
            self.logger.exception("validation failed")
            self.update_processing_status(dataset_version_id, DatasetStatusKey.ATAC, DatasetConversionStatus.FAILED)
            raise ValidationAtacFailed(errors=[str(e)]) from None

        if errors:
            # if the validation fails, update the processing status and raise a ValidationAtacFailed exception
            self.update_processing_status(dataset_version_id, DatasetStatusKey.ATAC, DatasetConversionStatus.FAILED)
            raise ValidationAtacFailed(errors=errors)

        # Changes to processing only happen during a migration. Only hash the files if the migration is set to true
        in_migration = os.environ.get("MIGRATION", "").lower() == "true"
        if in_migration:
            # check if the new fragment is the same as the old fragment
            fragment_unchanged = self.hash_file(local_fragment_filename) == self.hash_file(fragment_file)
        else:
            fragment_unchanged = False

        # fragment file to avoid uploading the same file multiple times
        # if the fragment file is unchanged from a migration or the fragment file is already ingested, use the old fragment.
        if fragment_unchanged or (self.business_logic.is_already_ingested(manifest.atac_fragment) and not in_migration):
            # get the artifact id of the old fragment, and add it to the new dataset
            artifact_name = str(manifest.atac_fragment).split("/")[-1]
            artifact = self.business_logic.database_provider.get_artifact_by_uri_suffix(artifact_name)
            self.business_logic.database_provider.add_artifact_to_dataset_version(dataset_version_id, artifact.id)
            # get the artifact id of the old fragment index, and add it to the new dataset
            artifact = self.business_logic.database_provider.get_artifact_by_uri_suffix(artifact_name + ".tbi")
            self.business_logic.database_provider.add_artifact_to_dataset_version(dataset_version_id, artifact.id)
            self.update_processing_status(dataset_version_id, DatasetStatusKey.ATAC, DatasetConversionStatus.COPIED)
        else:
            fragment_artifact_id = self.create_atac_artifact(
                fragment_file,
                DatasetArtifactType.ATAC_FRAGMENT,
                dataset_version_id,
                datasets_bucket,
            )
            self.create_atac_artifact(
                fragment_index_file,
                DatasetArtifactType.ATAC_INDEX,
                dataset_version_id,
                datasets_bucket,
                fragment_artifact_id,
            )
            self.update_processing_status(dataset_version_id, DatasetStatusKey.ATAC, DatasetConversionStatus.UPLOADED)
        self.logger.info("Processing completed successfully")
        return
