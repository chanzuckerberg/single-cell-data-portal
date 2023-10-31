import contextlib
import logging
import shutil

import requests

from backend.common.utils.corpora_constants import CorporaConstants
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
from backend.layers.thirdparty.uri_provider import UriProviderInterface

logger = logging.getLogger(__name__)


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

        # This is a bit ugly and should be done polymorphically instead, but Dropbox support will be dropped soon
        if file_url.scheme == "https":
            file_size = self.uri_provider.get_file_info(source_uri).size
            if file_size and file_size >= shutil.disk_usage("/")[2]:
                raise UploadFailed("Insufficient disk space.")
            download(file_url.url, local_path)
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


def download(
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
            logger.debug("Starting download.")
            for chunk in resp.iter_content(chunk_size=chunk_size):
                if chunk:
                    fp.write(chunk)
                    chunk_size = len(chunk)
                    logger.debug(f"chunk size: {chunk_size}")
