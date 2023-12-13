import contextlib
import logging

import requests

logger = logging.getLogger("processing")


def download(
    uri: str,
    local_path: str,
    chunk_size: int = 10 * 2**20,
) -> None:
    """
    Download a file from a URI and update the processing_status upload fields in the database

    :param uri: The URI of the file to be downloaded.
    :param local_path: The local name of the file be downloaded.
    :param chunk_size: The size of downloaded data to copy to memory before saving to disk.

    :return: The current dataset processing status.
    """

    with contextlib.suppress(Exception), requests.get(uri, stream=True) as resp:
        resp.raise_for_status()
        with open(local_path, "wb") as fp:
            logger.debug("Starting download.")
            for chunk in resp.iter_content(chunk_size=chunk_size):
                if chunk:
                    fp.write(chunk)
                    chunk_size = len(chunk)
                    logger.debug(f"chunk size: {chunk_size}")
