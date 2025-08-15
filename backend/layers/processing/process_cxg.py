from typing import Optional

from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetStatusKey,
    DatasetVersionId,
)
from backend.layers.processing.h5ad_data_file import H5ADDataFile
from backend.layers.processing.logger import logit
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.thirdparty.s3_provider import S3ProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface


class ProcessCxg(ProcessingLogic):
    """
    Base class for handling the `Process CXG` step of the step function.
    This will:
    1. Download the labeled h5ad artifact from S3 (uploaded by DownloadAndValidate)
    2. Convert to cxg
    3. Upload the cxg artifact (a directory) to S3
    If this step completes successfully, the handle_success lambda will be invoked
    If this step fails, the handle_failures lambda will be invoked
    """

    def __init__(
        self,
        business_logic: BusinessLogicInterface,
        uri_provider: UriProviderInterface,
        s3_provider: S3ProviderInterface,
    ) -> None:
        super().__init__()
        self.business_logic = business_logic
        self.uri_provider = uri_provider
        self.s3_provider = s3_provider

    def process(
        self,
        dataset_version_id: DatasetVersionId,
        artifact_bucket: str,
        cellxgene_bucket: str,
        is_reprocess=False,
    ):
        """
        1. Download the labeled dataset from the artifact bucket
        2. Convert the labeled dataset to CXG
        3. Upload the CXG to the cellxgene bucket
        :param dataset_version_id:
        :param artifact_bucket:
        :param cellxgene_bucket:
        :param is_reprocess: flag indicating whether this job is reprocessing an existing cxg in-place
        :return:
        """

        labeled_h5ad_filename = "local.h5ad"

        # Download the labeled dataset from the artifact bucket
        object_key = None
        current_artifacts = self.business_logic.get_dataset_artifacts(dataset_version_id)
        if is_reprocess:
            existing_h5ad = [artifact for artifact in current_artifacts if artifact.type == DatasetArtifactType.H5AD][0]
            if existing_h5ad:
                _, object_key = self.s3_provider.parse_s3_uri(existing_h5ad.uri)
        if object_key is None:
            key_prefix = self.get_key_prefix(dataset_version_id.id)
            object_key = f"{key_prefix}/{labeled_h5ad_filename}"
        self.download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

        # Get fragment file S3 URL
        fragment_file_path = None
        fragment_artifacts = [
            artifact for artifact in current_artifacts if artifact.type == DatasetArtifactType.ATAC_FRAGMENT
        ]
        if fragment_artifacts:
            fragment_file_path = fragment_artifacts[0].uri
            self.logger.info(f"Using fragment file S3 URL: {fragment_file_path}")

        # Convert the labeled dataset to CXG and upload it to the cellxgene bucket
        self.process_cxg(
            labeled_h5ad_filename,
            dataset_version_id,
            cellxgene_bucket,
            fragment_file_path,
            current_artifacts,
            is_reprocess,
        )

    @logit
    def make_cxg(self, local_filename, dataset_version_id: DatasetVersionId, fragment_file_path: Optional[str] = None):
        """
        Convert the uploaded H5AD file to the CXG format servicing the cellxgene Explorer.
        """

        cxg_output_container = local_filename.replace(".h5ad", ".cxg")
        try:
            h5ad_data_file = H5ADDataFile(local_filename)
            h5ad_data_file.to_cxg(
                cxg_output_container,
                sparse_threshold=25.0,
                dataset_version_id=dataset_version_id.id,
                fragment_artifact_id=fragment_file_path,
            )
        except Exception as ex:
            # TODO use a specialized exception
            msg = "CXG conversion failed."
            self.logger.exception(msg)
            raise RuntimeError(msg) from ex

        return cxg_output_container

    def delete_existing_cxg_files(self, s3_uri):
        """
        Delete all existing files in a given CXG S3 directory
        """
        bucket, prefix = self.s3_provider.parse_s3_uri(s3_uri)
        prefix = prefix.rstrip("/")
        self.logger.info(f"Deleting existing CXG files from s3://{bucket}/{prefix}/")
        self.s3_provider.delete_prefix(bucket, prefix)

    def copy_cxg_files_to_cxg_bucket(self, cxg_dir, s3_uri, clear_existing=False):
        """
        Copy cxg files to the cellxgene bucket (under the given object key) for access by the explorer
        """
        if clear_existing:
            self.delete_existing_cxg_files(s3_uri)
        self.s3_provider.upload_directory(cxg_dir, s3_uri)

    def process_cxg(
        self,
        local_filename,
        dataset_version_id,
        cellxgene_bucket,
        fragment_file_path=None,
        current_artifacts=None,
        is_reprocess=False,
    ):
        cxg_dir = self.convert_file(
            self.make_cxg, local_filename, dataset_version_id, fragment_file_path, DatasetStatusKey.CXG
        )
        s3_uri = None
        existing_cxg_artifacts = []
        if current_artifacts:
            existing_cxg_artifacts = [
                artifact for artifact in current_artifacts if artifact.type == DatasetArtifactType.CXG
            ]
            if existing_cxg_artifacts:
                s3_uri = existing_cxg_artifacts[0].uri

        if s3_uri is None:
            key_prefix = self.get_key_prefix(dataset_version_id.id)
            s3_uri = f"s3://{cellxgene_bucket}/{key_prefix}.cxg/"

        self.update_processing_status(dataset_version_id, DatasetStatusKey.CXG, DatasetConversionStatus.UPLOADING)
        self.copy_cxg_files_to_cxg_bucket(cxg_dir, s3_uri, clear_existing=is_reprocess)
        self.logger.info(f"Updating database with cxg artifact for dataset {dataset_version_id}. s3_uri is {s3_uri}")

        if not existing_cxg_artifacts:
            self.business_logic.add_dataset_artifact(dataset_version_id, DatasetArtifactType.CXG, s3_uri)

        self.update_processing_status(dataset_version_id, DatasetStatusKey.CXG, DatasetConversionStatus.UPLOADED)
