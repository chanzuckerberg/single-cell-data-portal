import json
import logging
import os
from dataclasses import dataclass, field
from typing import Dict, Optional

import pandas as pd
import tiledb
from pandas import DataFrame
from tiledb import Array

from backend.common.utils.s3_buckets import buckets
from backend.wmg.config import WmgConfig
from backend.wmg.data.tiledb import create_ctx

# Snapshot data artifact file/dir names
CELL_TYPE_ORDERINGS_FILENAME = "cell_type_orderings.json"
PRIMARY_FILTER_DIMENSIONS_FILENAME = "primary_filter_dimensions.json"
EXPRESSION_SUMMARY_CUBE_NAME = "expression_summary"
EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME = "expression_summary_default"
CELL_COUNTS_CUBE_NAME = "cell_counts"
MARKER_GENES_CUBE_NAME = "marker_genes"
DATASET_TO_GENE_IDS_FILENAME = "dataset_to_gene_ids.json"
FILTER_RELATIONSHIPS_FILENAME = "filter_relationships.json"
DATASET_METADATA_FILENAME = "dataset_metadata.json"

STACK_NAME = os.environ.get("REMOTE_DEV_PREFIX")

# root directory under which the data artifact exists
WMG_ROOT_DIR_PATH = STACK_NAME.strip("/") if STACK_NAME else ""

logger = logging.getLogger("wmg")

###################################### PUBLIC INTERFACE #################################


@dataclass
class WmgSnapshot:
    """
    All of the data artifacts the WMG API depends upon to perform its functions, versioned by "snapshot_identifier".
    These are read from data artifacts, per the relative file names, above.
    """

    snapshot_identifier: Optional[str] = field(default=None)

    # TileDB array containing expression summary statistics (expressed gene count, non-expressed mean,
    # etc.) aggregated by multiple cell metadata dimensions and genes. See the full schema at
    # backend/wmg/data/schemas/cube_schema.py.
    expression_summary_cube: Optional[Array] = field(default=None)

    # TileDB array containing expression summary statistics optimized for querying with no
    # secondary filters selected.
    # See the full schema at backend/wmg/data/schemas/cube_schema_default.py.
    expression_summary_default_cube: Optional[Array] = field(default=None)

    # TileDB array containing the precomputed marker genes.
    # See the full schema at backend/wmg/data/schemas/marker_gene_cube_schema.py.
    marker_genes_cube: Optional[Array] = field(default=None)

    # TileDB array containing the total cell counts (expressed gene count, non-expressed mean, etc.) aggregated by
    # multiple cell metadata dimensions (but no gene dimension). See the full schema at
    # backend/wmg/data/schemas/cube_schema.py.
    cell_counts_cube: Optional[Array] = field(default=None)

    # Pandas DataFrame containing per-tissue ordering of cell types.
    # Columns are "tissue_ontology_term_id", "cell_type_ontology_term_id", "order"
    cell_type_orderings: Optional[DataFrame] = field(default=None)

    # precomputed list of ids for all gene and tissue ontology term ids per organism
    primary_filter_dimensions: Optional[Dict] = field(default=None)

    # precomputed filter relationships graph
    filter_relationships: Optional[Dict] = field(default=None)

    # dataset metadata dictionary
    dataset_metadata: Optional[Dict] = field(default=None)


# Cached data
cached_snapshot: Optional[WmgSnapshot] = None


def load_snapshot(
    *,
    snapshot_schema_version: str,
    explicit_snapshot_id_to_load: Optional[str] = None,
    read_versioned_snapshot: bool = True,
) -> WmgSnapshot:
    """
    Loads and caches the snapshot identified by the snapshot schema version and a snapshot id.

    By default, this functions loads the latest snapshot id for a given schema version.
    If an explicit snapshot id is given, it will load that snapshot id. Note, that if an explicit
    snapshot id is given, it should be associated with the given schema version.

    The snapshot representation is cached in memory. Therefore, multiple calls to this function
    will simply return the cached snapshot if there isn't a newer snapshot id.
    """
    global cached_snapshot

    should_reload, snapshot_id = _should_reload_snapshot(
        snapshot_schema_version=snapshot_schema_version,
        explicit_snapshot_id_to_load=explicit_snapshot_id_to_load,
        read_versioned_snapshot=read_versioned_snapshot,
    )

    if should_reload:
        cached_snapshot = _load_snapshot(
            snapshot_schema_version=snapshot_schema_version,
            snapshot_id=snapshot_id,
            read_versioned_snapshot=read_versioned_snapshot,
        )
    # TODO: Once the pipeline generates the V2 snapshot, this can be removed
    if cached_snapshot.dataset_metadata is None:
        cached_snapshot.build_dataset_metadata_dict()
    return cached_snapshot


###################################### PRIVATE INTERFACE #################################
def _get_latest_snapshot_id(*, snapshot_schema_version: str, read_versioned_snapshot: bool) -> str:
    """
    Get latest snapshot id for a given snapshot schema version
    """
    data_schema_dir_path = _get_wmg_snapshot_schema_dir_path(
        snapshot_schema_version=snapshot_schema_version,
        read_versioned_snapshot=read_versioned_snapshot,
    )
    file_name = "latest_snapshot_identifier"

    key_path = f"{data_schema_dir_path}/{file_name}" if data_schema_dir_path else file_name

    latest_snapshot_id = _read_value_at_s3_key(key_path)
    return latest_snapshot_id


def _get_wmg_snapshot_schema_dir_path(*, snapshot_schema_version: str, read_versioned_snapshot: bool) -> str:
    """
    Get path to a particular snapshot schema version.
    """
    data_schema_dir_path = f"snapshots/{snapshot_schema_version}"

    if WMG_ROOT_DIR_PATH:
        data_schema_dir_path = f"{WMG_ROOT_DIR_PATH}/{data_schema_dir_path}"

    if read_versioned_snapshot:
        return data_schema_dir_path
    return WMG_ROOT_DIR_PATH


def _get_wmg_snapshot_dir_path(*, snapshot_schema_version: str, snapshot_id: str, read_versioned_snapshot: bool) -> str:
    """
    Get path to the snapshot id directory for a the given snapshot schema version.
    """
    data_schema_dir_path = _get_wmg_snapshot_schema_dir_path(
        snapshot_schema_version=snapshot_schema_version,
        read_versioned_snapshot=read_versioned_snapshot,
    )

    snapshot_id_dir_path = f"{data_schema_dir_path}/{snapshot_id}" if data_schema_dir_path else snapshot_id
    return snapshot_id_dir_path


def _get_wmg_snapshot_dir_s3_uri(snapshot_dir_path: str) -> str:
    """
    Get the s3 uri of the snapshot directory
    """
    wmg_config = WmgConfig()

    return os.path.join("s3://", wmg_config.bucket, snapshot_dir_path)


def _load_snapshot(*, snapshot_schema_version: str, snapshot_id: str, read_versioned_snapshot: bool) -> WmgSnapshot:
    snapshot_dir_path = _get_wmg_snapshot_dir_path(
        snapshot_schema_version=snapshot_schema_version,
        snapshot_id=snapshot_id,
        read_versioned_snapshot=read_versioned_snapshot,
    )

    cell_type_orderings = _load_cell_type_order(snapshot_dir_path)
    primary_filter_dimensions = _load_primary_filter_data(snapshot_dir_path)
    filter_relationships = _load_filter_graph_data(snapshot_dir_path)

    # TODO: Once the pipeline generates the V2 snapshot, the ternary can be removed
    dataset_metadata = _load_dataset_metadata(snapshot_dir_path) if snapshot_schema_version != "v1" else None

    snapshot_base_uri = _get_wmg_snapshot_dir_s3_uri(snapshot_dir_path)
    logger.info(f"Loading WMG snapshot from directory path {snapshot_dir_path} into URI {snapshot_base_uri}")

    # TODO: Okay to keep TileDB arrays open indefinitely? Is it faster than re-opening each request?
    #  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell
    #  -data-portal/2134
    return WmgSnapshot(
        snapshot_identifier=snapshot_id,
        expression_summary_cube=_open_cube(f"{snapshot_base_uri}/{EXPRESSION_SUMMARY_CUBE_NAME}"),
        expression_summary_default_cube=_open_cube(f"{snapshot_base_uri}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}"),
        marker_genes_cube=_open_cube(f"{snapshot_base_uri}/{MARKER_GENES_CUBE_NAME}"),
        cell_counts_cube=_open_cube(f"{snapshot_base_uri}/{CELL_COUNTS_CUBE_NAME}"),
        cell_type_orderings=cell_type_orderings,
        primary_filter_dimensions=primary_filter_dimensions,
        filter_relationships=filter_relationships,
        dataset_metadata=dataset_metadata,
    )


def _open_cube(cube_uri) -> Array:
    return tiledb.open(cube_uri, ctx=create_ctx(json.loads(WmgConfig().tiledb_config_overrides)))


def _load_cell_type_order(snapshot_dir_path: str) -> DataFrame:
    key_path = f"{snapshot_dir_path}/{CELL_TYPE_ORDERINGS_FILENAME}"
    return pd.read_json(_read_value_at_s3_key(key_path))


def _load_primary_filter_data(snapshot_dir_path: str) -> Dict:
    key_path = f"{snapshot_dir_path}/{PRIMARY_FILTER_DIMENSIONS_FILENAME}"
    return json.loads(_read_value_at_s3_key(key_path))


def _load_dataset_metadata(snapshot_dir_path: str) -> Dict:
    key_path = f"{snapshot_dir_path}/{DATASET_METADATA_FILENAME}"
    return json.loads(_read_value_at_s3_key(key_path))


def _load_filter_graph_data(snapshot_dir_path: str) -> str:
    try:
        key_path = f"{snapshot_dir_path}/{FILTER_RELATIONSHIPS_FILENAME}"
        return json.loads(_read_value_at_s3_key(key_path))
    except Exception:
        logger.warning(
            f"{_get_wmg_snapshot_dir_s3_uri(snapshot_dir_path)}/{FILTER_RELATIONSHIPS_FILENAME} could not be loaded"
        )
        return None


def _read_value_at_s3_key(key_path: str):
    """
    Read value at an s3 key
    """
    s3 = buckets.portal_resource

    wmg_config = WmgConfig()
    wmg_config.load()

    s3obj = s3.Object(wmg_config.bucket, key_path)
    return s3obj.get()["Body"].read().decode("utf-8").strip()


def _should_reload_snapshot(
    *,
    snapshot_schema_version: str,
    explicit_snapshot_id_to_load: Optional[str] = None,
    read_versioned_snapshot: bool,
) -> tuple[bool, str]:
    """
    Returns a pair: (<should_reload>, <snapshot_id>) where <should_reload> is a boolean indicating
    whether then in-memory snapshot should be reloaded and <snapshot_id> is the id of the snapshot that
    the in-memory data structure represents.
    """
    snapshot_id = explicit_snapshot_id_to_load or _get_latest_snapshot_id(
        snapshot_schema_version=snapshot_schema_version, read_versioned_snapshot=read_versioned_snapshot
    )

    if cached_snapshot is None:
        logger.info(f"Loading snapshot id: {snapshot_id}")
        return (True, snapshot_id)
    elif snapshot_id != cached_snapshot.snapshot_identifier:
        logger.info(
            f"Reloading snapshot. Detected snapshot id update from cached "
            f"snapshot id: {cached_snapshot.snapshot_identifier} to {snapshot_id}"
        )
        return (True, snapshot_id)
    else:
        logger.debug(f"No need to load snapshot. Using cached snapshot id: {cached_snapshot.snapshot_identifier}")
        return (False, snapshot_id)
