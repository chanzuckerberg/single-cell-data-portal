import logging
import os
import uuid
from typing import Tuple

import boto3
import tiledb
from tiledb import Array

from backend.wmg.config import WmgConfig
from backend.wmg.data.tiledb import fast_config, create_ctx

logger = logging.getLogger('wmg')

# Cached cube
cube = None
latest_snapshot_identifier_s3obj = None
latest_snapshot_identifier = None


def load_cube() -> Tuple[Array, str]:
    """
    Loads and caches the WMG cube (TileDB Array). Reloads the cube if the latest_snapshot_identifier S3 object has
    been updated.
    @return: TileDB Array and latest snapshot identifier, as a Tuple
    """

    global cube, latest_snapshot_identifier

    if _update_latest_snapshot_identifier() or cube is None:
        # TODO: Okay to keep open indefinitely? Is it faster than re-opening each request?
        cube_uri = build_cube_uri(WmgConfig().bucket, latest_snapshot_identifier)
        logger.info(f"Opening WMG cube at {cube_uri}")
        cube = _open_cube(cube_uri)

    return cube, latest_snapshot_identifier


def _open_cube(cube_uri) -> Array:
    return tiledb.open(cube_uri, ctx=create_ctx(fast_config()))


# TODO: Worth doing this on a thread, continuously, rather than on-demand, in order to proactively open a new cube (
#  and maybe warm the TileDB cache?) before a user needs to query it
def _update_latest_snapshot_identifier() -> bool:
    global latest_snapshot_identifier_s3obj, latest_snapshot_identifier

    # if latest_snapshot_identifier_s3obj is None:
        # s3 = boto3.resource('s3')
        # latest_snapshot_identifier_s3obj = s3.Object(WmgConfig().bucket, 'latest_snapshot_identifier')

    # latest_snapshot_identifier_object.reload() # necessary?
    # new_snapshot_identifier = latest_snapshot_identifier_s3obj.get()['Body'].read().decode('utf-8')
    new_snapshot_identifier = latest_snapshot_identifier or uuid.uuid4().hex

    if new_snapshot_identifier != latest_snapshot_identifier:
        logger.info(f'detected snapshot update from {latest_snapshot_identifier} to {new_snapshot_identifier}')
        latest_snapshot_identifier = new_snapshot_identifier
        return True
    else:
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f'latest snapshot identifier={latest_snapshot_identifier}')
        return False


def build_latest_snapshot_uri(data_root_uri: str):
    return os.path.join(data_root_uri, 'latest_snapshot_identifier')


def build_cube_uri(bucket: str, snapshot_identifier: str):
    return os.path.join("s3://", bucket, snapshot_identifier, 'cube')



