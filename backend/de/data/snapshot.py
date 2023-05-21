import json
import logging
import os
from dataclasses import dataclass
from typing import Dict, Optional

import tiledb
from tiledb import Array

from backend.common.utils.s3_buckets import buckets
from backend.de.config import DeConfig
from backend.de.data.tiledb import create_ctx

# Snapshot data artifact file/dir names
EXPRESSION_SUMMARY_CUBE_NAMES = [
    "expression_summary__default",
    "expression_summary__dataset_id",
    "expression_summary__disease_ontology_term_id",
    "expression_summary__sex_ontology_term_id",
    "expression_summary__self_reported_ethnicity_ontology_term_id",
]
CELL_COUNTS_CUBE_NAME = "cell_counts"
FILTER_RELATIONSHIPS_FILENAME = f"filter_relationships.json"
CARDINALITY_PER_DIMENSION_FILENAME = f"cardinality_per_dimension.json"

logger = logging.getLogger("de")


@dataclass
class DeSnapshot:
    """
    All of the data artifacts the DE API depends upon to perform its functions, versioned by "snapshot_identifier".
    These are read from data artifacts, per the relative file names, above.
    """

    snapshot_identifier: str

    expression_summary_cubes: Dict[str,Array]

    cell_counts_cube: Array

    # precomputed filter relationships graph
    filter_relationships: Dict

    cardinality_per_dimension: Dict

    def __hash__(self):
        return hash(None)  # hash is not used for DeSnapshot

# Cached data
cached_snapshot: Optional[DeSnapshot] = None


def load_snapshot() -> DeSnapshot:
    global cached_snapshot
    if new_snapshot_identifier := _update_latest_snapshot_identifier():
        cached_snapshot = _load_snapshot(new_snapshot_identifier)
        cached_snapshot.build_dataset_metadata_dict()
    return cached_snapshot


def _load_snapshot(new_snapshot_identifier) -> DeSnapshot:
    filter_relationships = _load_filter_graph_data(new_snapshot_identifier)
    cardinality_per_dimension = _load_cardinality_per_dimension_data(new_snapshot_identifier)

    snapshot_base_uri = _build_snapshot_base_uri(new_snapshot_identifier)
    logger.info(f"Loading DE snapshot at {snapshot_base_uri}")

    expression_summary_cubes = {}
    for name in EXPRESSION_SUMMARY_CUBE_NAMES:
        dim = name.split("__")[-1]
        expression_summary_cubes[dim] = _open_cube(f"{snapshot_base_uri}/{name}")

    return DeSnapshot(
        snapshot_identifier=new_snapshot_identifier,
        expression_summary_cubes=expression_summary_cubes,
        cell_counts_cube=_open_cube(f"{snapshot_base_uri}/{CELL_COUNTS_CUBE_NAME}"),
        filter_relationships=filter_relationships,
        cardinality_per_dimension=cardinality_per_dimension,
    )


def _open_cube(cube_uri) -> Array:
    return tiledb.open(cube_uri, ctx=create_ctx(json.loads(DeConfig().tiledb_config_overrides)))


def _load_filter_graph_data(snapshot_identifier: str):
    try:
        return json.loads(_read_s3obj(f"{snapshot_identifier}/{FILTER_RELATIONSHIPS_FILENAME}"))
    except Exception:
        return None

def _load_cardinality_per_dimension_data(snapshot_identifier: str):
    try:
        return json.loads(_read_s3obj(f"{snapshot_identifier}/{CARDINALITY_PER_DIMENSION_FILENAME}"))
    except Exception:
        return None
    
def _read_s3obj(relative_path: str) -> str:
    s3 = buckets.portal_resource
    de_config = DeConfig()
    de_config.load()
    prefixed_relative_path = os.path.join(_build_data_path_prefix(), relative_path or "")
    s3obj = s3.Object(DeConfig().bucket, prefixed_relative_path)
    return s3obj.get()["Body"].read().decode("utf-8").strip()

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


def _build_data_path_prefix():
    de_config = DeConfig()
    rdev_prefix = os.environ.get("REMOTE_DEV_PREFIX")
    if rdev_prefix:
        return rdev_prefix.strip("/")
    elif "data_path_prefix" in de_config.config:
        return de_config.data_path_prefix
    else:
        return ""


def _build_snapshot_base_uri(snapshot_identifier: str):
    de_config = DeConfig()
    data_path_prefix = _build_data_path_prefix()
    return os.path.join(
        "s3://",
        de_config.bucket,
        data_path_prefix,
        snapshot_identifier,
    )
