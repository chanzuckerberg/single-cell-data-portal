import logging
import threading

import requests

from ..common.corpora_orm import DbDatasetProcessingStatus
from ..common.entities import Dataset
from ..common.utils.db_utils import db_session_manager
from ..common.utils.math_utils import MB

logger = logging.getLogger(__name__)


class ProgressTracker:
    def __init__(self, file_size: int):
        self.file_size: int = file_size
        self._progress: int = 0
        self.lock: threading.Lock = threading.Lock()
        self.complete: threading.Event = threading.Event()

    def progress(self):
        with self.lock:
            return self._progress / self.file_size

    def update(self, progress):
        with self.lock:
            self._progress += progress


def uploader(url: str, local_path: str, tracker: ProgressTracker, chunk_size: int):
    """

    :param url: The URL of the file to be uploaded.
    :param local_path: The local name of the file be uploaded
    :param tracker: Tracks information about the progress of the upload.
    :return:
    """
    try:
        with requests.get(url, stream=True) as resp:
            resp.raise_for_status()
            with open(local_path, "wb") as fp:
                for chunk in resp.iter_content(chunk_size=chunk_size):
                    if chunk:
                        fp.write(chunk)
                        chunk_size = len(chunk)
                        tracker.update(chunk_size)
                        logger.debug(f"chunk size: {chunk_size}")
            logger.info("Upload Complete!")
    finally:
        tracker.complete.set()


def update_progress(upload_uuid: str, tracker: ProgressTracker, frequency: float = 3):
    """
    Update the progress of an upload to the database using the tracker.

    :param upload_uuid: The uuid of the upload_progress row.
    :param tracker: Tracks information about the progress of the upload.
    :param frequency: The frequency in which the database is updated
    :return:
    """

    def _update():
        progress = tracker.progress()
        with db_session_manager() as manager:
            manager.session.query(DbDatasetProcessingStatus).filter(DbDatasetProcessingStatus.id == upload_uuid).update(
                {DbDatasetProcessingStatus.upload_progress: progress}
            )
            manager.session.commit()
            logger.info(f"progress: {100 * progress:.2f} %")

    while not tracker.complete.wait(frequency):
        _update()
    _update()  # Make sure the progress is update once the upload is complete


def upload(dataset_uuid: str, url: str, local_path: str, file_size: int, chunk_size: int = 10 * MB, update_frequency=3):
    """

    :param dataset_uuid: The uuid of the dataset the upload will be associated with.
    :param url: The URL of the file to be uploaded.
    :param local_path: The local name of the file be uploaded
    :param file_size: The size of the file in bytes.
    :param chunk_size: The size of downloaded data to copy to memory before saving to disk.
    :return:
    """
    with db_session_manager():
        dataset = Dataset.get(dataset_uuid)
        processing_status = dataset.new_processing_status()
        dataset.update(processing_status=processing_status)
        upload_uuid = dataset.processing_status.id
    progress_tracker = ProgressTracker(file_size)

    progress_thread = threading.Thread(
        target=update_progress,
        kwargs=dict(upload_uuid=upload_uuid, tracker=progress_tracker, frequency=update_frequency),
    )
    progress_thread.start()

    upload_thread = threading.Thread(
        target=uploader, kwargs=dict(url=url, local_path=local_path, tracker=progress_tracker, chunk_size=chunk_size)
    )
    upload_thread.start()
    upload_thread.join()
    progress_thread.join()
