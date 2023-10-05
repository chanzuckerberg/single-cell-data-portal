import logging
import os
import subprocess

from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetStatusKey,
    DatasetVersionId,
)
from backend.layers.processing.logging_utils import logit
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.thirdparty.s3_provider import S3ProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface

logger: logging.Logger = logging.getLogger("processing")


class ProcessSeurat(ProcessingLogic):
    """
    Base class for handling the `Process Seurat` step of the step function.
    This will:
    1. Determine if a Seurat conversion is possible (this is set by the previous step - DownloadAndValidate)
    2. Download the labeled h5ad artifact from S3 (uploaded by DownloadAndValidate)
    3. Convert to RDS
    4. Upload the RDS artifact to S3
    If this step completes successfully, and ProcessCxg is completed, the handle_success lambda will be invoked
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

    def process(self, dataset_id: DatasetVersionId, artifact_bucket: str, datasets_bucket: str):
        """
        1. Download the labeled dataset from the artifact bucket
        2. Convert it to Seurat format
        3. Upload the Seurat file to the artifact bucket
        :param dataset_id:
        :param artifact_bucket:
        :param datasets_bucket:
        :return:
        """

        # If the validator previously marked the dataset as rds_status.SKIPPED, do not start the Seurat processing
        dataset = self.business_logic.get_dataset_version(dataset_id)

        if dataset is None:
            raise Exception("Dataset not found")  # TODO: maybe improve

        if dataset.status.rds_status == DatasetConversionStatus.SKIPPED:
            self.logger.info("Skipping Seurat conversion")
            return

        labeled_h5ad_filename = "local.h5ad"

        key_prefix = self.get_key_prefix(dataset_id.id)
        object_key = f"{key_prefix}/{labeled_h5ad_filename}"
        self.download_from_s3(artifact_bucket, object_key, labeled_h5ad_filename)

        seurat_filename = self.convert_file(
            self.make_seurat,
            labeled_h5ad_filename,
            "Failed to convert dataset to Seurat format.",
            dataset_id,
            DatasetStatusKey.RDS,
        )

        self.create_artifact(
            seurat_filename,
            DatasetArtifactType.RDS,
            key_prefix,
            dataset_id,
            artifact_bucket,
            DatasetStatusKey.RDS,
            datasets_bucket=datasets_bucket,
        )

    @logit
    def make_seurat(self, local_filename):
        """
        Create a Seurat rds file from the AnnData file.
        """
        try:
            completed_process = subprocess.run(
                ["Rscript", "-e", "\"installed.packages()[, c('Package', 'Version')]\""], capture_output=True
            )
            logger.debug({"stdout": completed_process.stdout, "args": completed_process.args})

            subprocess.run(
                [
                    "Rscript",
                    os.path.join(os.path.abspath(os.path.dirname(__file__)), "make_seurat.R"),
                    local_filename,
                ],
                capture_output=True,
                check=True,
            )
        except subprocess.CalledProcessError as ex:
            msg = f"Seurat conversion failed: {ex.output} {ex.stderr}"
            self.logger.exception(msg)
            raise RuntimeError(msg) from ex

        return local_filename.replace(".h5ad", ".rds")
