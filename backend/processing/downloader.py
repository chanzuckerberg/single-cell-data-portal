import contextlib
import logging
import shutil
import threading

import requests

from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import DatasetStatusKey, DatasetUploadStatus, DatasetVersionId
from backend.processing.exceptions import UploadFailed

logger = logging.getLogger(__name__)


class ProgressTracker:
    def __init__(self, file_size: int):
        self.file_size: int = file_size
        self._progress: int = 0
        self.progress_lock: threading.Lock = threading.Lock()  # prevent concurrent access of ProgressTracker._progress
        self.stop_updater: threading.Event = threading.Event()  # Stops the update_progress thread
        self.stop_downloader: threading.Event = threading.Event()  # Stops the downloader threads
        self.error: Exception = None  # Track errors
        self.tombstoned: bool = False  # Track if dataset tombstoned

    def progress(self):
        with self.progress_lock:
            return self._progress / self.file_size

    def update(self, progress):
        with self.progress_lock:
            self._progress += progress

    def cancel(self):
        self.tombstoned = True
        self.stop_downloader.set()
        self.stop_updater.set()


class NoOpProgressTracker:
    """
    This progress tracker should be used if file_size isn't available.
    It will return a progress of 1 if no errors occurred
    during the download (i.e. if self.error was never set), otherwise it will return 0.
    """

    def __init__(self) -> None:
        self.progress_lock: threading.Lock = threading.Lock()  # prevent concurrent access of ProgressTracker._progress
        self.stop_updater: threading.Event = threading.Event()  # Stops the update_progress thread
        self.stop_downloader: threading.Event = threading.Event()  # Stops the downloader threads
        self.error: Exception = None  # Track errors
        self.tombstoned: bool = False  # Track if dataset tombstoned

    def progress(self):
        if self.error:
            return 0
        else:
            return 1

    def update(self, progress):
        pass

    def cancel(self):
        self.tombstoned = True
        self.stop_downloader.set()
        self.stop_updater.set()


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

        # if file_size is not None:
        #     progress_tracker = ProgressTracker(file_size)
        # else:
        #     progress_tracker = NoOpProgressTracker()

        with contextlib.suppress(Exception):
            self.download_file(url, local_path, chunk_size)
            # TODO: maybe add a check on the file size

        self.business_logic.update_dataset_version_status(
            dataset_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADED
        )

        # progress_thread.join()  # Wait for the progress thread to complete
        # if progress_tracker.tombstoned:
        #     raise ProcessingCancelled
        # if progress_tracker.error:
        #     status = {
        #         "upload_status": UploadStatus.FAILED,
        #         "upload_message": str(progress_tracker.error),
        #     }
        #     _processing_status_updater(processing_status, status)

        # status_dict = processing_status.to_dict()
        # if processing_status.upload_status == UploadStatus.FAILED:
        #     logger.error(f"Upload failed: {status_dict}")
        #     raise ProcessingFailed(status_dict)
        # else:
        #     return status_dict
