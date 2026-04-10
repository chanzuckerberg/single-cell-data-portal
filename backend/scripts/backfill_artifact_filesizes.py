"""
Backfill filesize for existing DatasetArtifact rows.

Fetches the size of each artifact from S3 using concurrent HEAD requests and
writes the result back to the database. Safe to re-run: rows with a non-NULL
filesize are skipped.

Usage (from repo root, with DB and AWS credentials available):
    python -m backend.scripts.backfill_artifact_filesizes

Optional env vars:
    BACKFILL_WORKERS   number of concurrent S3 threads (default: 50)
    BACKFILL_DRY_RUN   if set to "1", print counts but make no DB writes
"""

import logging
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

from backend.layers.persistence.orm import DatasetArtifactTable
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.thirdparty.s3_provider import S3Provider

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

MAX_WORKERS = int(os.environ.get("BACKFILL_WORKERS", 50))
DRY_RUN = os.environ.get("BACKFILL_DRY_RUN", "") == "1"


def _fetch_filesize(s3: S3Provider, artifact_id: str, uri: str):
    """Return (artifact_id, filesize) or (artifact_id, None) on failure."""
    try:
        size = s3.get_file_size(uri)
        return artifact_id, size
    except Exception as exc:
        logger.warning("S3 HEAD failed for %s (%s): %s", artifact_id, uri, exc)
        return artifact_id, None


def backfill():
    db = DatabaseProvider()
    s3 = S3Provider()

    with db._manage_session() as session:
        rows = (
            session.query(DatasetArtifactTable.id, DatasetArtifactTable.uri)
            .filter(DatasetArtifactTable.filesize.is_(None))
            .all()
        )

    total = len(rows)
    logger.info("Found %d artifact(s) with NULL filesize%s", total, " (DRY RUN)" if DRY_RUN else "")
    if total == 0 or DRY_RUN:
        return

    succeeded = 0
    failed = 0

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(_fetch_filesize, s3, str(row.id), row.uri): str(row.id) for row in rows}
        for future in as_completed(futures):
            artifact_id, filesize = future.result()
            if filesize is None:
                failed += 1
                continue
            with db._manage_session() as session:
                session.query(DatasetArtifactTable).filter_by(id=artifact_id).update(
                    {"filesize": filesize}, synchronize_session=False
                )
            succeeded += 1
            if succeeded % 100 == 0:
                logger.info("Progress: %d/%d updated, %d failed", succeeded, total, failed)

    logger.info("Done. %d updated, %d failed out of %d total.", succeeded, failed, total)
    if failed:
        logger.warning("%d artifacts could not be resolved — re-run to retry.", failed)


if __name__ == "__main__":
    backfill()
