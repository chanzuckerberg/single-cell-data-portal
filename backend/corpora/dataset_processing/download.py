import logging
import shutil
import threading

import requests
from sqlalchemy import inspect

from backend.corpora.common.corpora_orm import DbDatasetProcessingStatus, UploadStatus, ProcessingStatus
from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.utils.math_utils import MB
from backend.corpora.dataset_processing.exceptions import ProcessingFailed, ProcessingCancelled

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


def downloader(url: str, local_path: str, tracker: ProgressTracker, chunk_size: int):
    """
    Download the file pointed at by the URL to the local path.

    :param url: The URL of the file to be downloaded.
    :param local_path: The local name of the file to be downloaded
    :param tracker: Tracks information about the progress of the download.
    :param chunk_size: The size of downloaded data to copy to memory before saving to disk.
    :return:
    """
    try:
        with requests.get(url, stream=True) as resp:
            resp.raise_for_status()
            with open(local_path, "wb") as fp:
                logger.info("Starting download.")
                for chunk in resp.iter_content(chunk_size=chunk_size):
                    if tracker.stop_downloader.is_set():
                        logger.warning("Download ended early!")
                        return
                    elif chunk:
                        fp.write(chunk)
                        chunk_size = len(chunk)
                        tracker.update(chunk_size)
                        logger.debug(f"chunk size: {chunk_size}")
    except Exception as ex:
        tracker.error = ex
        logger.error(f"Download Failed for {url}")
    finally:
        tracker.stop_updater.set()


def updater(processing_status: DbDatasetProcessingStatus, tracker: ProgressTracker, frequency: float):
    """
    Update the progress of an upload to the database using the tracker.

    :param processing_status: the SQLAlchemy ProcessingStatus object
    :param tracker: Tracks information about the progress of the upload.
    :param frequency: The frequency in which the database is updated in seconds
    :return:
    """
    db_session = inspect(processing_status).session

    def _update():
        progress = tracker.progress()
        dataset = Dataset.get(db_session, processing_status.dataset_id, include_tombstones=True)
        if dataset.tombstone:
            tracker.cancel()
            return

        elif tracker.stop_updater.is_set():
            if progress > 1:
                tracker.stop_downloader.set()
                message = "The expected file size is smaller than the downloaded file size."
                status = {
                    "upload_progress": progress,
                    "upload_message": message,
                    "upload_status": UploadStatus.FAILED,
                }
            elif progress < 1:
                message = "The expected file size is greater than the actual file size."
                status = {
                    "upload_progress": progress,
                    "upload_message": message,
                    "upload_status": UploadStatus.FAILED,
                }
            elif progress == 1:
                status = {
                    "upload_progress": progress,
                    "upload_status": UploadStatus.UPLOADED,
                }
        else:
            status = {"upload_progress": progress}
        logger.debug("Updating processing_status")
        _processing_status_updater(processing_status, status)
        db_session.commit()

    try:
        while not tracker.stop_updater.wait(frequency):
            _update()
        _update()  # Make sure the progress is updated once the download is complete
    finally:
        tracker.stop_downloader.set()


def _processing_status_updater(processing_status: DbDatasetProcessingStatus, updates: dict):
    for key, value in updates.items():
        setattr(processing_status, key, value)


def download(
    dataset_uuid: str,
    url: str,
    local_path: str,
    file_size: int,
    chunk_size: int = 10 * MB,
    update_frequency=3,
) -> dict:
    """
    Download a file from a url and update the processing_status upload fields in the database

    :param dataset_uuid: The uuid of the dataset the download will be associated with.
    :param url: The URL of the file to be downloaded.
    :param local_path: The local name of the file be downloaded.
    :param file_size: The size of the file in bytes.
    :param chunk_size: Forwarded to downloader thread
    :param update_frequency: The frequency in which to update the database in seconds.

    :return: The current dataset processing status.
    """
    with db_session_manager() as session:
        logger.info("Setting up download.")
        logger.info(f"file_size: {file_size}")
        if file_size >= shutil.disk_usage("/")[2]:
            raise ProcessingFailed("Insufficient disk space.")
        processing_status = Dataset.get(session, dataset_uuid).processing_status
        processing_status.upload_status = UploadStatus.UPLOADING
        processing_status.upload_progress = 0
        progress_tracker = ProgressTracker(file_size)

        progress_thread = threading.Thread(
            target=updater,
            kwargs=dict(processing_status=processing_status, tracker=progress_tracker, frequency=update_frequency),
        )
        progress_thread.start()
        download_thread = threading.Thread(
            target=downloader,
            kwargs=dict(url=url, local_path=local_path, tracker=progress_tracker, chunk_size=chunk_size),
        )
        download_thread.start()
        download_thread.join()  # Wait for the download thread to complete
        progress_thread.join()  # Wait for the progress thread to complete
        if progress_tracker.tombstoned:
            raise ProcessingCancelled
        if progress_tracker.error:
            status = {
                "upload_status": UploadStatus.FAILED,
                "upload_message": str(progress_tracker.error),
                "processing_status": ProcessingStatus.FAILURE,
            }
            _processing_status_updater(processing_status, status)

        status_dict = processing_status.to_dict()
        if processing_status.upload_status == UploadStatus.FAILED:
            logger.error(f"Upload failed: {status_dict}")
            raise ProcessingFailed(status_dict)
        else:
            return status_dict
