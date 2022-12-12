#!/usr/bin/env python3
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetStatusKey,
    DatasetVersionId,
)
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.thirdparty.s3_provider import S3ProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface
from backend.portal.pipeline.processing.h5ad_data_file import H5ADDataFile


class ProcessCxg(ProcessingLogic):
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

    def process(self, dataset_id: DatasetVersionId, artifact_bucket: str, cellxgene_bucket: str):
        """
        1. Download the labeled dataset from the artifact bucket
        2. Convert the labeled dataset to CXG
        3. Upload the CXG to the cellxgene bucket
        :param dataset_id:
        :param artifact_bucket:
        :param cellxgene_bucket:
        :return:
        """

        labeled_h5ad_filename = "local.h5ad"

        # Download the labeled dataset from the artifact bucket
        bucket_prefix = self.get_bucket_prefix(dataset_id.id)
        object_key = f"{bucket_prefix}/{labeled_h5ad_filename}"
        self.download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

        # Convert the labeled dataset to CXG and upload it to the cellxgene bucket
        self.process_cxg(labeled_h5ad_filename, dataset_id, cellxgene_bucket)

    def make_cxg(self, local_filename):
        """
        Convert the uploaded H5AD file to the CXG format servicing the cellxgene Explorer.
        """

        cxg_output_container = local_filename.replace(".h5ad", ".cxg")
        try:
            h5ad_data_file = H5ADDataFile(local_filename, var_index_column_name="feature_name")
            h5ad_data_file.to_cxg(cxg_output_container, sparse_threshold=25.0)
        except Exception as ex:
            # TODO use a specialized exception
            msg = "CXG conversion failed."
            self.logger.exception(msg)
            raise RuntimeError(msg) from ex

        return cxg_output_container

    def copy_cxg_files_to_cxg_bucket(self, cxg_dir, s3_uri):
        """
        Copy cxg files to the cellxgene bucket (under the given object key) for access by the explorer
        """
        self.s3_provider.upload_directory(cxg_dir, s3_uri)

    def process_cxg(self, local_filename, dataset_id, cellxgene_bucket):
        cxg_dir = self.convert_file(
            self.make_cxg, local_filename, "Issue creating cxg.", dataset_id, DatasetStatusKey.CXG
        )
        bucket_prefix = self.get_bucket_prefix(dataset_id.id)
        s3_uri = f"s3://{cellxgene_bucket}/{bucket_prefix}.cxg/"
        self.update_processing_status(dataset_id, DatasetStatusKey.CXG, DatasetConversionStatus.UPLOADING)
        self.copy_cxg_files_to_cxg_bucket(cxg_dir, s3_uri)
        # TODO: this is where we set explorer_url, but we might not need it anymore
        # metadata = {
        #     "explorer_url": join(
        #         DEPLOYMENT_STAGE_TO_URL[os.environ["DEPLOYMENT_STAGE"]],
        #         dataset_id + ".cxg",
        #         "",
        #     )
        # }
        self.business_logic.add_dataset_artifact(dataset_id, DatasetArtifactType.CXG, s3_uri)
        self.logger.info(f"Updating database with cxg artifact for dataset {dataset_id}. s3_uri is {s3_uri}")
        self.update_processing_status(dataset_id, DatasetStatusKey.CXG, DatasetConversionStatus.UPLOADED)
