import logging
import threading

import requests

from ..common.corpora_orm import DbDatasetProcessingStatus
from ..common.entities import Dataset
from ..common.utils.db_utils import db_session_manager

logger = logging.getLogger(__name__)


class ProgressTracker:
    def __init__(self, total: int):
        self.total: int = total
        self._progress: int = 0
        self.lock: threading.Lock = threading.Lock()
        self.complete: threading.Event = threading.Event()

    def progress(self):
        with self.lock:
            _progress = self._progress
        return _progress / self.total

    def update(self, progress):
        with self.lock:
            self._progress += progress


class Upload(threading.Thread):
    def __init__(self, url: str, local_path: str, tracker: ProgressTracker):
        threading.Thread.__init__(self)
        self.tracker = tracker
        self.url = url
        self.local_path = local_path

    def run(self):
        tracker = self.tracker
        with requests.get(self.url, stream=True) as resp:
            resp.raise_for_status()
            with open(self.local_path, "wb") as fp:
                for chunk in resp.iter_content():
                    if chunk:
                        fp.write(chunk)
                        chunk_size = len(chunk)
                        tracker.update(chunk_size)
                        logger.debug(f"chunk size: {chunk_size}")
        tracker.complete.set()


class Progress(threading.Thread):
    def __init__(self, dataset_uuid, tracker: ProgressTracker):
        threading.Thread.__init__(self)
        self.tracker: ProgressTracker = tracker
        with db_session_manager():
            self.uuid = Dataset.get(dataset_uuid).processing_status.id

    def run(self):
        tracker = self.tracker
        while not tracker.complete.wait(3):
            progress = tracker.progress()
            with db_session_manager() as session:
                session.query(DbDatasetProcessingStatus).filter(DbDatasetProcessingStatus.id == self.uuid).update(
                    {DbDatasetProcessingStatus.upload_progress: progress}
                )
                session.commit()
                logger.debug(f"update {progress}")
