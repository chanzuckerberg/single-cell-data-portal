import contextlib
import logging
import shutil

import requests

from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import DatasetStatusKey, DatasetUploadStatus, DatasetVersionId
from backend.layers.processing.exceptions import UploadFailed

logger = logging.getLogger("processing")


class Downloader:

    business_logic: BusinessLogicInterface

    def __init__(self, business_logic: BusinessLogicInterface) -> None:
        self.business_logic = business_logic

    def download_file(self, url: str, local_path: str, chunk_size: int):
        """
        Download the file pointed at by the URL to the local path.

        :param url: The URL of the file to be downloaded.
        :param local_path: The local name of the file to be downloaded
        :param chunk_size: The size of downloaded data to copy to memory before saving to disk.
        :return:
        """
        with requests.get(url, stream=True) as resp:
            resp.raise_for_status()
            with open(local_path, "wb") as fp:
                logger.info("Starting download.")
                for chunk in resp.iter_content(chunk_size=chunk_size):
                    if chunk:
                        fp.write(chunk)
                        chunk_size = len(chunk)
                        logger.debug(f"chunk size: {chunk_size}")

    def download(
        self,
        dataset_id: DatasetVersionId,
        url: str,
        local_path: str,
        file_size: int,
        chunk_size: int = 10 * 2**20,
        update_frequency=3,
    ) -> None:
        """
        Download a file from a url and update the processing_status upload fields in the database

        :param dataset_id: The uuid of the dataset the download will be associated with.
        :param url: The URL of the file to be downloaded.
        :param local_path: The local name of the file be downloaded.
        :param file_size: The size of the file in bytes.
        :param chunk_size: Forwarded to downloader thread
        :param update_frequency: The frequency in which to update the database in seconds.

        :return: The current dataset processing status.
        """
        logger.info("Setting up download.")
        logger.info(f"file_size: {file_size}")

        if file_size and file_size >= shutil.disk_usage("/")[2]:
            raise UploadFailed("Insufficient disk space.")

        self.business_logic.update_dataset_version_status(
            dataset_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING
        )
        # TODO: set upload_progress to 0

        with contextlib.suppress(Exception):
            self.download_file(url, local_path, chunk_size)
            # TODO: maybe add a check on the file size

        self.business_logic.update_dataset_version_status(
            dataset_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADED
        )
