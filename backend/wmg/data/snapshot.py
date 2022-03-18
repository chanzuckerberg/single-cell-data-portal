import logging
import os
from collections import namedtuple
from dataclasses import dataclass
from typing import Optional

import pandas as pd
import tiledb
from pandas import DataFrame
from tiledb import Array

from backend.corpora.common.utils.s3_buckets import buckets
from backend.wmg.config import WmgConfig
from backend.wmg.data.tiledb import create_ctx

logger = logging.getLogger("wmg")

@dataclass
class WmgSnapshot:
    """
    All of the data artifacts the WMG API depends upon to perform its functions, versioned by "snapshot_identifier".
    """
    snapshot_identifier: str
    expression_summary_cube: Array
    cell_counts_cube: Array
    cell_type_orderings: DataFrame


# Cached data
cached_snapshot: Optional[WmgSnapshot] = None


def load_snapshot() -> WmgSnapshot:
    """
    Loads and caches the WMG snapshot. Reloads the snapshot data if the latest_snapshot_identifier S3 object has
    been updated.
    @return: WmgSnapshot object
    """

    global cached_snapshot

    if new_snapshot_identifier := _update_latest_snapshot_identifier():
        cached_snapshot = _load_snapshot(new_snapshot_identifier)
    return cached_snapshot


def _load_snapshot(new_snapshot_identifier) -> WmgSnapshot:
    snapshot_base_uri = _build_snapshot_base_uri(WmgConfig().bucket, new_snapshot_identifier)
    logger.info(f"Loading WMG snapshot at {snapshot_base_uri}")
    # TODO: Okay to keep TileDB arrays open indefinitely? Is it faster than re-opening each request?
    #  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell
    #  -data-portal/2134
    return WmgSnapshot(snapshot_identifier=new_snapshot_identifier,
                       expression_summary_cube=_open_cube(f"{snapshot_base_uri}/expression_summary"),
                       cell_counts_cube=_open_cube(f"{snapshot_base_uri}/cell_counts"),
                       cell_type_orderings=_load_cell_type_order())


def _open_cube(cube_uri) -> Array:
    return tiledb.open(cube_uri, ctx=create_ctx(tiledb_mem_gb=float(WmgConfig().tiledb_mem_gb)))


def _load_cell_type_order() -> DataFrame:
    return pd.read_json(_read_s3obj("cell_type_orderings.json"))


def _read_s3obj(relative_path: str) -> str:
    s3 = buckets.portal_resource
    s3obj = s3.Object(WmgConfig().bucket, relative_path)
    return s3obj.get()["Body"].read().decode("utf-8").strip()

# TODO: Worth doing this on a thread, continuously, rather than on-demand, in order to proactively load an updated
#  snapshot (and maybe warm the TileDB caches?) before a user needs to query the data:
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2134
def _update_latest_snapshot_identifier() -> Optional[str]:
    global cached_snapshot

    new_snapshot_identifier = _read_s3obj("latest_snapshot_identifier")

    if cached_snapshot is None:
        logger.info(f"using latest snapshot {new_snapshot_identifier}")
        return new_snapshot_identifier
    elif new_snapshot_identifier != cached_snapshot.snapshot_identifier:
        logger.info(f"detected snapshot update from {cached_snapshot.snapshot_identifier} to {new_snapshot_identifier}")
        return new_snapshot_identifier
    else:
        logger.debug(f"latest snapshot identifier={cached_snapshot.snapshot_identifier}")
        return None


def _build_snapshot_base_uri(bucket: str, snapshot_identifier: str):
    return os.path.join("s3://", bucket, snapshot_identifier)
