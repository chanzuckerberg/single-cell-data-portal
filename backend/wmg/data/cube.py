import logging
import os
from collections import namedtuple
from typing import Tuple, Optional

import tiledb
from tiledb import Array

from backend.corpora.common.utils.s3_buckets import buckets
from backend.wmg.config import WmgConfig
from backend.wmg.data.tiledb import create_ctx

logger = logging.getLogger("wmg")

WmgCubes = namedtuple("cubes", ['expression_summary_cube', 'cell_counts_cube', 'snapshot_identifier'])

# Cached cubes
cubes: WmgCubes = None
latest_snapshot_identifier_s3obj = None


def load_cubes() -> WmgCubes:
    """
    Loads and caches the WMG cube (TileDB Array). Reloads the cube if the latest_snapshot_identifier S3 object has
    been updated.
    @return: TileDB Array and latest snapshot identifier, as a Tuple
    """

    global cubes

    if new_snapshot_identifier := _update_latest_snapshot_identifier():
        # TODO: Okay to keep open indefinitely? Is it faster than re-opening each request?
        #  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell
        #  -data-portal/2134
        snapshot_base_uri = build_snapshot_base_uri(WmgConfig().bucket, new_snapshot_identifier)
        logger.info(f"Opening WMG cube at {snapshot_base_uri}")
        cubes = _open_cubes(snapshot_base_uri, new_snapshot_identifier)
    return cubes


def _open_cubes(snapshot_base_uri, new_snapshot_identifier) -> WmgCubes:
    expression_summary = _open_cube(f"{snapshot_base_uri}/expression_summary")
    cell_counts = _open_cube(f"{snapshot_base_uri}/cell_counts")
    return WmgCubes(expression_summary, cell_counts, new_snapshot_identifier)


def _open_cube(cube_uri) -> Array:
    return tiledb.open(cube_uri, ctx=create_ctx(tiledb_mem_gb=float(WmgConfig().tiledb_mem_gb)))


# TODO: Worth doing this on a thread, continuously, rather than on-demand, in order to proactively open a new cube (
#  and maybe warm the TileDB cache?) before a user needs to query it:
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2134
def _update_latest_snapshot_identifier() -> Optional[str]:
    global latest_snapshot_identifier_s3obj, cubes

    if latest_snapshot_identifier_s3obj is None:
        s3 = buckets.portal_resource
        latest_snapshot_identifier_s3obj = s3.Object(WmgConfig().bucket, "latest_snapshot_identifier")

    latest_snapshot_identifier_s3obj.reload()  # necessary?
    new_snapshot_identifier = latest_snapshot_identifier_s3obj.get()["Body"].read().decode("utf-8").strip()

    if cubes is None:
        logger.info(f"using latest snapshot {new_snapshot_identifier}")
        return new_snapshot_identifier
    elif new_snapshot_identifier != cubes.snapshot_identifier:
        logger.info(f"detected snapshot update from {cubes.snapshot_identifier} to {new_snapshot_identifier}")
        return new_snapshot_identifier
    else:
        logger.debug(f"latest snapshot identifier={cubes.snapshot_identifier}")
        return None


def build_latest_snapshot_uri(data_root_uri: str):
    return os.path.join(data_root_uri, "latest_snapshot_identifier")


def build_snapshot_base_uri(bucket: str, snapshot_identifier: str):
    return os.path.join("s3://", bucket, snapshot_identifier)
