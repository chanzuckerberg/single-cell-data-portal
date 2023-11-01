import json
import os
from typing import Any, Dict

import scanpy

from backend.common.utils.corpora_constants import CorporaConstants
from backend.common.utils.dl_sources.uri import DownloadFailed
from backend.common.utils.math_utils import MB
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    DatasetArtifactType,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetVersionId,
)
from backend.layers.processing.exceptions import UploadFailed
from backend.layers.processing.logger import logit
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.thirdparty.s3_provider_interface import S3ProviderInterface
from backend.layers.thirdparty.step_function_provider import StepFunctionProvider
from backend.layers.thirdparty.uri_provider import UriProviderInterface


class ProcessDownload(ProcessingLogic):
    """
    Base class for handling the `Download and Validate` step of the step function.
    This will:
    1. Download the original artifact from the provided URI
    2. estimate memory requirements
    4. Upload a copy of the original artifact (raw.h5ad)

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

    @logit
    def download_from_source_uri(self, source_uri: str, local_path: str) -> str:
        """Given a source URI, download it to local_path.
        Handles fixing the url so it downloads directly.
        """
        file_url = self.uri_provider.parse(source_uri)
        if not file_url:
            raise ValueError(f"Malformed source URI: {source_uri}")
        try:
            file_url.download(local_path)
        except DownloadFailed as e:
            raise UploadFailed(f"Failed to download file from source URI: {source_uri}") from e
        return local_path

    # TODO: after upgrading to Python 3.9, replace this with removeprefix()
    @staticmethod
    def remove_prefix(string: str, prefix: str) -> str:
        if string.startswith(prefix):
            return string[len(prefix) :]
        else:
            return string

    def create_batch_job_definition_parameters(self, local_filename: str, dataset_id) -> Dict[str, Any]:
        MEMORY_MODIFIER = 1.1  # add 10% overhead
        MEMORY_PER_VCPU = 4000
        MIN_MEMORY_MB = 4000
        MIN_VCPU = 1
        # The largest machine we are allocating is r5ad.2xlarge. This machine has 64GB of memory and 16 vCPUs.
        MAX_MEMORY_MB = 64000
        MAX_VCPU = 16
        SWAP_MEMORY_MB = 300000

        adata = scanpy.read_h5ad(local_filename)

        # Note: this is a rough estimate of the uncompressed size of the dataset. This method avoid loading the entire
        # dataset into memory.
        uncompressed_size_MB = adata.n_obs * adata.n_vars // MB
        estimated_memory_MB = max([uncompressed_size_MB * MEMORY_MODIFIER, MIN_MEMORY_MB])
        if estimated_memory_MB > MAX_MEMORY_MB:
            estimated_memory_MB = MAX_MEMORY_MB
            estimated_vcpus = MAX_VCPU
            max_swap = SWAP_MEMORY_MB
        else:
            estimated_vcpus = max(estimated_memory_MB // MEMORY_PER_VCPU, MIN_VCPU)
            max_swap = 0
        stack_name = os.environ["REMOTE_DEV_PREFIX"].replace("/", "")
        job_definition_name = f"dp-{os.environ['DEPLOYMENT_STAGE']}-{stack_name}-ingest-process-{dataset_id}"
        return {  # Using PascalCase to match the Batch API
            "JobDefinitionName": job_definition_name,
            "Vcpus": estimated_vcpus,
            "Memory": estimated_memory_MB,
            "LinuxParameters": {
                "Swappiness": 60,
                "MaxSwap": max_swap,
            },
        }

    def process(self, dataset_id: DatasetVersionId, dataset_uri: str, artifact_bucket: str, sfn_task_token: str):
        """
        1. Download the original dataset
        2. Upload the labeled dataset to the artifact bucket
        :param dataset_id:
        :param dataset_uri:
        :param artifact_bucket:
        :param sfn_task_token: use to report back the memory requirements
        :return:
        """

        self.update_processing_status(dataset_id, DatasetStatusKey.PROCESSING, DatasetProcessingStatus.PENDING)

        # Download the original dataset from Dropbox
        local_filename = self.download_from_source_uri(
            source_uri=dataset_uri,
            local_path=CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME,
        )

        response = self.create_batch_job_definition_parameters(local_filename, dataset_id.id)
        self.logger.info(response)

        key_prefix = self.get_key_prefix(dataset_id.id)
        # Upload the original dataset to the artifact bucket
        self.update_processing_status(dataset_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)
        self.create_artifact(
            local_filename,
            DatasetArtifactType.RAW_H5AD,
            key_prefix,
            dataset_id,
            artifact_bucket,
            DatasetStatusKey.H5AD,
        )
        self.update_processing_status(dataset_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADED)

        sfn_client = StepFunctionProvider().client
        sfn_client.send_task_success(taskToken=sfn_task_token, output=json.dumps(response))
