import contextlib
import json
import os
import shutil
from math import ceil
from typing import Any, Dict

import requests
import scanpy

from backend.common.utils.corpora_constants import CorporaConstants
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

MEMORY_MODIFIER = 1.1  # add 10% overhead
MEMORY_PER_VCPU = 4000
MIN_MEMORY_MB = 4000
MIN_VCPU = 1
# The largest machine we are allocating is r5ad.2xlarge. This machine has 64GB of memory and 16 vCPUs.
MAX_MEMORY_MB = 64000
MAX_VCPU = 16
SWAP_MEMORY_MB = 300000


class ProcessDownload(ProcessingLogic):
    """
    Base class for handling the `Download` step of the step function.
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

        # This is a bit ugly and should be done polymorphically instead, but Dropbox support will be dropped soon
        if file_url.scheme == "https":
            file_size = self.uri_provider.get_file_info(source_uri).size
            if file_size and file_size >= shutil.disk_usage("/")[2]:
                raise UploadFailed("Insufficient disk space.")
            self.download(file_url.url, local_path)
        elif file_url.scheme == "s3":
            bucket_name = file_url.netloc
            key = self.remove_prefix(file_url.path, "/")
            self.download_from_s3(
                bucket_name=bucket_name,
                object_key=key,
                local_filename=local_path,
            )
        else:
            raise ValueError(f"Download for URI scheme '{file_url.scheme}' not implemented")
        return local_path

    # TODO: after upgrading to Python 3.9, replace this with removeprefix()
    @staticmethod
    def remove_prefix(string: str, prefix: str) -> str:
        if string.startswith(prefix):
            return string[len(prefix) :]
        else:
            return string

    @staticmethod
    def get_job_definion_name(dataset_version_id: str) -> str:
        if os.getenv("REMOTE_DEV_PREFIX"):
            stack_name = os.environ["REMOTE_DEV_PREFIX"].replace("/", "")
            prefix = f"{os.environ['DEPLOYMENT_STAGE']}-{stack_name}"
        else:
            prefix = f"{os.environ['DEPLOYMENT_STAGE']}"
        job_definition_name = f"dp-{prefix}-ingest-process-{dataset_version_id}"
        return job_definition_name

    @staticmethod
    def estimate_resource_requirements(
        adata: scanpy.AnnData,
        memory_modifier: float = MEMORY_MODIFIER,
        min_memory_MB: int = MIN_MEMORY_MB,
        min_vcpu: int = MIN_VCPU,
        max_memory_MB: int = MAX_MEMORY_MB,
        max_vcpu: int = MAX_VCPU,
        swap_memory_MB: int = SWAP_MEMORY_MB,
        memory_per_vcpu: int = MEMORY_PER_VCPU,
    ) -> Dict[str, int]:
        """
        Estimate the resource requirements for a given dataset

        :param adata: The datasets AnnData object
        :param memory_modifier: A multiplier to increase/decrease the memory requirements by
        :param min_memory_MB: The minimum amount of memory to allocate.
        :param min_vcpu: The minimum number of vCPUs to allocate.
        :param max_memory_MB: The maximum amount of memory to allocate.
        :param max_vcpu: The maximum number of vCPUs to allocate.
        :param memory_per_vcpu: The amount of memory to allocate per vCPU.
        :param swap_memory_MB:
        :return: A dictionary containing the resource requirements
        """
        # Note: this is a rough estimate of the uncompressed size of the dataset. This method avoid loading the entire
        # dataset into memory.
        uncompressed_size_MB = adata.n_obs * adata.n_vars / MB
        estimated_memory_MB = max([int(ceil(uncompressed_size_MB * memory_modifier)), min_memory_MB])
        if estimated_memory_MB > max_memory_MB:
            estimated_memory_MB = max_memory_MB
            estimated_vcpus = max_vcpu
            max_swap = swap_memory_MB
        else:
            estimated_vcpus = max([int(ceil(estimated_memory_MB / memory_per_vcpu)), min_vcpu])
            max_swap = 0

        return {"Vcpus": estimated_vcpus, "Memory": estimated_memory_MB, "MaxSwap": max_swap}

    def create_batch_job_definition_parameters(self, local_filename: str, dataset_version_id: str) -> Dict[str, Any]:
        adata = scanpy.read_h5ad(local_filename, backed="r")
        batch_resources = self.estimate_resource_requirements(adata)
        job_definition_name = self.get_job_definion_name(dataset_version_id)

        return {  # Using PascalCase to match the Batch API
            "JobDefinitionName": job_definition_name,
            "Vcpus": batch_resources["Vcpus"],
            "Memory": batch_resources["Memory"],
            "LinuxParameters": {
                "Swappiness": 60,
                "MaxSwap": batch_resources["MaxSwap"],
            },
        }

    def process(
        self, dataset_version_id: DatasetVersionId, dataset_uri: str, artifact_bucket: str, sfn_task_token: str
    ):
        """
        Process the download step of the step function

        :param dataset_version_id:
        :param dataset_uri:
        :param artifact_bucket:
        :param sfn_task_token: use to report back the memory requirements
        :return:
        """

        self.update_processing_status(dataset_version_id, DatasetStatusKey.PROCESSING, DatasetProcessingStatus.PENDING)

        # Download the original dataset from Dropbox
        local_filename = self.download_from_source_uri(
            source_uri=dataset_uri,
            local_path=CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME,
        )

        response = self.create_batch_job_definition_parameters(local_filename, dataset_version_id.id)
        self.logger.info(response)

        key_prefix = self.get_key_prefix(dataset_version_id.id)
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

        sfn_client = StepFunctionProvider().client
        sfn_client.send_task_success(taskToken=sfn_task_token, output=json.dumps(response))

    def download(
        self,
        url: str,
        local_path: str,
        chunk_size: int = 10 * 2**20,
    ) -> None:
        """
        Download a file from a url and update the processing_status upload fields in the database

        :param url: The URL of the file to be downloaded.
        :param local_path: The local name of the file be downloaded.
        :param chunk_size: The size of downloaded data to copy to memory before saving to disk.

        :return: The current dataset processing status.
        """

        with contextlib.suppress(Exception), requests.get(url, stream=True) as resp:
            resp.raise_for_status()
            with open(local_path, "wb") as fp:
                self.logger.debug("Starting download.")
                for chunk in resp.iter_content(chunk_size=chunk_size):
                    if chunk:
                        fp.write(chunk)
                        chunk_size = len(chunk)
                        self.logger.debug(f"chunk size: {chunk_size}")
