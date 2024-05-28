import json
import logging
import os
from dataclasses import dataclass, field
from typing import Dict, Optional

import pandas as pd
import tiledb
from ddtrace import tracer
from pandas import DataFrame
from tiledb import Array

from backend.common.census_cube.config import CensusCubeConfig
from backend.common.census_cube.data.constants import CENSUS_CUBE_SNAPSHOT_FS_CACHE_ROOT_PATH
from backend.common.census_cube.data.tiledb import create_ctx
from backend.common.utils.s3_buckets import buckets

# Snapshot data artifact file/dir names
CELL_TYPE_ORDERINGS_FILENAME = "cell_type_orderings.json"
PRIMARY_FILTER_DIMENSIONS_FILENAME = "primary_filter_dimensions.json"
EXPRESSION_SUMMARY_CUBE_NAME = "expression_summary"
EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME = "expression_summary_default"
CELL_COUNTS_CUBE_NAME = "cell_counts"
EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAME = "expression_summary_diffexp"
EXPRESSION_SUMMARY_DIFFEXP_SIMPLE_CUBE_NAME = "expression_summary_diffexp_simple"
CELL_COUNTS_DIFFEXP_CUBE_NAME = "cell_counts_diffexp"
MARKER_GENES_CUBE_NAME = "marker_genes"
FILTER_RELATIONSHIPS_FILENAME = "filter_relationships.json"
DATASET_METADATA_FILENAME = "dataset_metadata.json"
CELL_TYPE_ANCESTORS_FILENAME = "cell_type_ancestors.json"

STACK_NAME = os.environ.get("REMOTE_DEV_PREFIX")

# root directory under which the data artifact exists
CENSUS_CUBE_ROOT_DIR_PATH = STACK_NAME.strip("/") if STACK_NAME else ""

DEPLOYMENT_STAGE = os.environ.get("DEPLOYMENT_STAGE", "")
SNAPSHOT_FS_ROOT_PATH = CENSUS_CUBE_SNAPSHOT_FS_CACHE_ROOT_PATH if (DEPLOYMENT_STAGE != "test") else None

logger = logging.getLogger("wmg")

###################################### PUBLIC INTERFACE #################################


@dataclass
class CensusCubeSnapshot:
    """
    All of the data artifacts the WMG API depends upon to perform its functions, versioned by "snapshot_identifier".
    These are read from data artifacts, per the relative file names, above.
    """

    snapshot_identifier: Optional[str] = field(default=None)

    # The version of the snapshot schema that is loaded.
    snapshot_schema_version: Optional[str] = field(default=None)

    # TileDB array containing expression summary statistics (expressed gene count, non-expressed mean,
    # etc.) aggregated by multiple cell metadata dimensions and genes. See the full schema at
    # backend/wmg/data/schemas/cube_schema.py.

    expression_summary_cube: Optional[Array] = field(default=None)

    # TileDB array containing expression summary statistics optimized for querying with no
    # secondary filters selected.
    # See the full schema at backend/wmg/data/schemas/cube_schema_default.py.
    expression_summary_default_cube: Optional[Array] = field(default=None)

    # TileDB arrays containing the precomputed expression summary statistics needed for
    # differential expression.
    diffexp_expression_summary_cubes: Optional[Dict[str, Array]] = field(default=None)

    # TileDB array containing the precomputed marker genes.
    # See the full schema at backend/wmg/data/schemas/marker_gene_cube_schema.py.
    marker_genes_cube: Optional[Array] = field(default=None)

    # TileDB array containing the total cell counts (expressed gene count, non-expressed mean, etc.) aggregated by
    # multiple cell metadata dimensions (but no gene dimension). See the full schema at
    # backend/wmg/data/schemas/cube_schema.py.
    # TODO: remove this in favor of cell_counts_df
    cell_counts_cube: Optional[Array] = field(default=None)

    # Dictionary of (cell type, tissue) tuples as keys and order as values.
    cell_type_orderings: Optional[dict] = field(default=None)

    # precomputed list of ids for all gene and tissue ontology term ids per organism
    primary_filter_dimensions: Optional[Dict] = field(default=None)

    # precomputed filter relationships graph
    filter_relationships: Optional[Dict] = field(default=None)

    # dataset metadata dictionary
    dataset_metadata: Optional[Dict] = field(default=None)

    # cell type ancestors pandas Series
    cell_type_ancestors: Optional[pd.Series] = field(default=None)

    # cell counts dataframe
    cell_counts_df: Optional[DataFrame] = field(default=None)

    # cell counts diffexp dataframe
    cell_counts_diffexp_df: Optional[DataFrame] = field(default=None)

    # expression summary diffexp cube
    expression_summary_diffexp_cube: Optional[Array] = field(default=None)

    # expression summary diffexp simple cube
    expression_summary_diffexp_simple_cube: Optional[Array] = field(default=None)


# Cached data
cached_snapshot: Optional[CensusCubeSnapshot] = None


@tracer.wrap(name="load_snapshot", service="wmg-api", resource="query", span_type="wmg-api")
def load_snapshot(
    *,
    snapshot_schema_version: str,
    explicit_snapshot_id_to_load: Optional[str] = None,
    snapshot_fs_root_path: Optional[str] = SNAPSHOT_FS_ROOT_PATH,
) -> CensusCubeSnapshot:
    """
    Loads and caches the snapshot identified by the snapshot schema version and a snapshot id.

    By default, this functions loads the latest snapshot id for a given schema version.
    If an explicit snapshot id is given, it will load that snapshot id. Note, that if an explicit
    snapshot id is given, it should be associated with the given schema version.

    The snapshot representation is cached in memory. Therefore, multiple calls to this function
    will simply return the cached snapshot if there isn't a newer snapshot id.

    Args:
        snapshot_schema_version (str): The version of the snapshot schema.
        explicit_snapshot_id_to_load (str, optional): The explicit snapshot id to load. Defaults to None.
        snapshot_fs_root_path (str, optional): The root path of the snapshot on the local filesystem. Defaults to None.

    Returns:
        CensusCubeSnapshot: The loaded snapshot.

    """
    global cached_snapshot

    if not (
        snapshot_fs_root_path
        and _local_disk_snapshot_is_valid(
            snapshot_fs_root_path=snapshot_fs_root_path,
            snapshot_schema_version=snapshot_schema_version,
            explicit_snapshot_id_to_load=explicit_snapshot_id_to_load,
        )
    ):
        snapshot_fs_root_path = None

    should_reload, snapshot_id = _should_reload_snapshot(
        snapshot_schema_version=snapshot_schema_version,
        explicit_snapshot_id_to_load=explicit_snapshot_id_to_load,
        snapshot_fs_root_path=snapshot_fs_root_path,
    )

    if should_reload:
        cached_snapshot = _load_snapshot(
            snapshot_schema_version=snapshot_schema_version,
            snapshot_id=snapshot_id,
            snapshot_fs_root_path=snapshot_fs_root_path,
        )

    return cached_snapshot


###################################### PRIVATE INTERFACE #################################
def _get_latest_snapshot_identifier_file_rel_path(snapshot_schema_version: str) -> str:
    """
    Get the relative path to the latest snapshot identifier file for a given snapshot schema version.

    Args:
        snapshot_schema_version (str): The version of the snapshot schema.

    Returns:
        str: The relative path to the latest snapshot identifier file.
    """
    data_schema_dir_path = _get_wmg_snapshot_schema_dir_rel_path(snapshot_schema_version)
    file_name = "latest_snapshot_identifier"

    return f"{data_schema_dir_path}/{file_name}"


def _get_latest_snapshot_id(snapshot_schema_version: str, snapshot_fs_root_path: Optional[str] = None) -> str:
    """
    Get latest snapshot id for a given snapshot schema version

    Args:
        snapshot_schema_version (str): The version of the snapshot schema.
        snapshot_fs_root_path (str, optional): The root path of the snapshot on the local filesystem. Defaults to None.

    Returns:
        str: The latest snapshot id for the given snapshot schema version.
    """
    rel_path = _get_latest_snapshot_identifier_file_rel_path(snapshot_schema_version)

    latest_snapshot_id = _read_wmg_data_file(rel_path, snapshot_fs_root_path)
    return latest_snapshot_id


def _get_wmg_snapshot_schema_dir_rel_path(snapshot_schema_version: str) -> str:
    """
    Get relative path to a particular snapshot schema version.

    That is, the relative path is the path that does not include root path where the data
    resides.

    Examples:

    1. An S3 fullpath to a snapshot schema version maybe s3://env-rdev-wmg/pr-6447/snapshots/v3.
    Here, "s3://env-rdev-wmg" is the S3 bucket. The S3 bucket is considered the "root" of the
    fullpath. Therefore, "pr-6447/snapshots/v3" would be the relative path returned by this function.

    2. A filesystem fullpath to a snapshot schema version maybe:
    /single-cell-data-portal/census_cube_snapshot_cache/snapshots/v3.
    Here, "/single-cell-data-portal/census_cube_snapshot_cache" is considered the "root" of the
    fullpath. Therefore "snapshots/v3" would be the relative path returned by
    this function.

    Args:
        snapshot_schema_version (str): The version of the snapshot schema.

    Returns:
        str: The relative path to the snapshot schema version.


    """
    data_schema_dir_rel_path = f"snapshots/{snapshot_schema_version}"

    if CENSUS_CUBE_ROOT_DIR_PATH:
        data_schema_dir_rel_path = f"{CENSUS_CUBE_ROOT_DIR_PATH}/{data_schema_dir_rel_path}"

    return data_schema_dir_rel_path


def _get_wmg_snapshot_rel_path(snapshot_schema_version: str, snapshot_id: str) -> str:
    """
    Get relative path to the snapshot id directory for a the given snapshot schema version.

    That is, the relative path is the path that does not include root path where the data
    resides.

    Examples:

    1. An S3 fullpath to a particular snapshot_id maybe: s3://env-rdev-wmg/pr-6447/snapshots/v3/1704754452.
    Here, "s3://env-rdev-wmg" is the S3 bucket. The S3 bucket is considered the "root" of the
    fullpath. Therefore, "pr-6447/snapshots/v3/1704754452" would be the relative path
    returned by this function.

    2. A filesystem fullpath to a particular snapshot_id maybe:
    /single-cell-data-portal/census_cube_snapshot_cache/snapshots/v3/1704754452.
    Here, "/single-cell-data-portal/census_cube_snapshot_cache" is considered the "root" of the
    fullpath. Therefore "snapshots/v3/1704754452" would be the relative path returned by
    this function.

    Args:
        snapshot_schema_version (str): The version of the snapshot schema.
        snapshot_id (str): The unique identifier of the snapshot.

    Returns:
        str: The relative path to the snapshot id directory for the given snapshot schema version.

    """
    data_schema_dir_rel_path = _get_wmg_snapshot_schema_dir_rel_path(snapshot_schema_version)

    snapshot_id_dir_rel_path = f"{data_schema_dir_rel_path}/{snapshot_id}"
    return snapshot_id_dir_rel_path


def _get_wmg_snapshot_fullpath(snapshot_rel_path: str, snapshot_fs_root_path: Optional[str] = None) -> str:
    """
    Return the full path of the snapshot on local disk or S3 URI of the snapshot.

    Examples:

    1. For snapshot on S3, this maybe: s3://env-rdev-wmg/pr-6447/snapshots/v3/1704754452.

    2. For snapshot on local disk, this maybe:
    /single-cell-data-portal/census_cube_snapshot_cache/snapshots/v3/1704754452

    Args:
        snapshot_rel_path (str): The relative path of the snapshot.
        snapshot_fs_root_path (Optional[str]): The root path of the snapshot in the filesystem. Defaults to None.

    Returns:
        str: The full path of the snapshot on local disk or S3 URI of the snapshot.
    """

    if snapshot_fs_root_path:
        return os.path.join(snapshot_fs_root_path, snapshot_rel_path)

    wmg_config = CensusCubeConfig()
    return os.path.join("s3://", wmg_config.bucket, snapshot_rel_path)


def _load_snapshot(
    *, snapshot_schema_version: str, snapshot_id: str, snapshot_fs_root_path: Optional[str] = None
) -> CensusCubeSnapshot:
    """
    Load a snapshot given its schema version, id, and root path in the filesystem.

    Args:
        snapshot_schema_version (str): The version of the snapshot schema.
        snapshot_id (str): The unique identifier of the snapshot.
        snapshot_fs_root_path (Optional[str]): The root path of the snapshot in the filesystem. Defaults to None.

    Returns:
        CensusCubeSnapshot: The loaded snapshot.
    """

    snapshot_rel_path = _get_wmg_snapshot_rel_path(snapshot_schema_version, snapshot_id)

    cell_type_orderings = _load_cell_type_order(snapshot_rel_path, snapshot_fs_root_path)
    primary_filter_dimensions = _load_primary_filter_data(snapshot_rel_path, snapshot_fs_root_path)
    filter_relationships = _load_filter_graph_data(snapshot_rel_path, snapshot_fs_root_path)
    cell_type_ancestors = _load_cell_type_ancestors(snapshot_rel_path, snapshot_fs_root_path)
    dataset_metadata = _load_dataset_metadata(snapshot_rel_path, snapshot_fs_root_path)

    snapshot_uri = _get_wmg_snapshot_fullpath(snapshot_rel_path, snapshot_fs_root_path)
    logger.info(f"Loading WMG snapshot from absolute path: {snapshot_uri}")

    # TODO: Okay to keep TileDB arrays open indefinitely? Is it faster than re-opening each request?
    #  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell
    #  -data-portal/2134
    cell_counts_cube = _open_cube(f"{snapshot_uri}/{CELL_COUNTS_CUBE_NAME}")
    cell_counts_diffexp_cube = _open_cube(f"{snapshot_uri}/{CELL_COUNTS_DIFFEXP_CUBE_NAME}")
    return CensusCubeSnapshot(
        snapshot_identifier=snapshot_id,
        snapshot_schema_version=snapshot_schema_version,
        expression_summary_cube=_open_cube(f"{snapshot_uri}/{EXPRESSION_SUMMARY_CUBE_NAME}"),
        expression_summary_default_cube=_open_cube(f"{snapshot_uri}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}"),
        marker_genes_cube=_open_cube(f"{snapshot_uri}/{MARKER_GENES_CUBE_NAME}"),
        cell_counts_cube=cell_counts_cube,
        cell_type_orderings=cell_type_orderings.set_index(["tissue_ontology_term_id", "cell_type_ontology_term_id"])[
            "order"
        ].to_dict(),
        primary_filter_dimensions=primary_filter_dimensions,
        filter_relationships=filter_relationships,
        dataset_metadata=dataset_metadata,
        cell_type_ancestors=pd.Series(cell_type_ancestors),
        cell_counts_df=cell_counts_cube.df[:],
        cell_counts_diffexp_df=cell_counts_diffexp_cube.df[:],
        expression_summary_diffexp_cube=_open_cube(f"{snapshot_uri}/{EXPRESSION_SUMMARY_DIFFEXP_CUBE_NAME}"),
        expression_summary_diffexp_simple_cube=_open_cube(
            f"{snapshot_uri}/{EXPRESSION_SUMMARY_DIFFEXP_SIMPLE_CUBE_NAME}"
        ),
    )


def _local_disk_snapshot_is_valid(
    *,
    snapshot_fs_root_path: str,
    snapshot_schema_version: str,
    explicit_snapshot_id_to_load: Optional[str] = None,
) -> bool:
    """
    Checks that the path on local disk contains valid WMG snapshot.

    Args:
        snapshot_fs_root_path (str): The root path of the snapshot on the local filesystem.
        snapshot_schema_version (str): The version of the snapshot schema.
        explicit_snapshot_id_to_load (Optional[str]): The explicit snapshot id to load. Defaults to None.

    Returns:
        bool: True if the path on local disk contains valid WMG snapshot, False otherwise.
    """
    latest_snapshot_identifier_file_rel_path = _get_latest_snapshot_identifier_file_rel_path(snapshot_schema_version)
    latest_snapshot_identifier_file_fullpath = os.path.join(
        snapshot_fs_root_path, latest_snapshot_identifier_file_rel_path
    )

    if not os.path.exists(latest_snapshot_identifier_file_fullpath):
        logger.warning(
            f"{latest_snapshot_identifier_file_fullpath} does not exist. Falling back to S3 to load WMG snapshot..."
        )
        return False
    else:
        if explicit_snapshot_id_to_load:
            snapshot_id = explicit_snapshot_id_to_load
        else:
            snapshot_id = _get_latest_snapshot_id(snapshot_schema_version, snapshot_fs_root_path)

        snapshot_rel_path = _get_wmg_snapshot_rel_path(snapshot_schema_version, snapshot_id)
        snapshot_full_path = _get_wmg_snapshot_fullpath(snapshot_rel_path, snapshot_fs_root_path)

        if not os.path.exists(snapshot_full_path):
            logger.warning(f"{snapshot_full_path} does not exist. Falling back to S3 to load WMG snapshot...")
            return False

    return True


def _open_cube(cube_uri) -> Array:
    return tiledb.open(cube_uri, ctx=create_ctx(json.loads(CensusCubeConfig().tiledb_config_overrides)))


def _load_cell_type_order(snapshot_rel_path: str, snapshot_fs_root_path: Optional[str] = None) -> DataFrame:
    rel_path = f"{snapshot_rel_path}/{CELL_TYPE_ORDERINGS_FILENAME}"
    return pd.read_json(_read_wmg_data_file(rel_path, snapshot_fs_root_path))


def _load_primary_filter_data(snapshot_rel_path: str, snapshot_fs_root_path: Optional[str] = None) -> Dict:
    rel_path = f"{snapshot_rel_path}/{PRIMARY_FILTER_DIMENSIONS_FILENAME}"
    return json.loads(_read_wmg_data_file(rel_path, snapshot_fs_root_path))


def _load_dataset_metadata(snapshot_rel_path: str, snapshot_fs_root_path: Optional[str] = None) -> Dict:
    rel_path = f"{snapshot_rel_path}/{DATASET_METADATA_FILENAME}"
    return json.loads(_read_wmg_data_file(rel_path, snapshot_fs_root_path))


def _load_cell_type_ancestors(snapshot_rel_path: str, snapshot_fs_root_path: Optional[str] = None) -> Dict:
    rel_path = f"{snapshot_rel_path}/{CELL_TYPE_ANCESTORS_FILENAME}"
    return json.loads(_read_wmg_data_file(rel_path, snapshot_fs_root_path))


def _load_filter_graph_data(snapshot_rel_path: str, snapshot_fs_root_path: Optional[str] = None) -> str:
    try:
        rel_path = f"{snapshot_rel_path}/{FILTER_RELATIONSHIPS_FILENAME}"
        return json.loads(_read_wmg_data_file(rel_path, snapshot_fs_root_path))
    except Exception:
        logger.warning(
            f"{_get_wmg_snapshot_fullpath(snapshot_rel_path)}/{FILTER_RELATIONSHIPS_FILENAME} could not be loaded"
        )
        return None


def _read_wmg_data_file(rel_path: str, snapshot_fs_root_path: Optional[str] = None) -> str:
    """
    Read file from local disk if snapshot_fs_root_path is provided. Otherwise, read from S3.

    When reading from S3, the 'rel_path' argument is the S3 key of the object to read.
    When reading from the local filesystem, 'rel_path' is suffixed to the 'snapshot_fs_root_path'
    to derive a fullpath of the file to read.

    Args:
        rel_path (str): The relative path of the file to read.
        snapshot_fs_root_path (Optional[str]): The root path of the snapshot in the filesystem. Defaults to None.

    Returns:
        str: The content of the file as a string.
    """
    if snapshot_fs_root_path:
        full_path = os.path.join(snapshot_fs_root_path, rel_path)
        with open(full_path, encoding="utf-8") as f:
            data = f.read().strip()
        return data

    return _read_value_at_s3_key(key_path=rel_path)


def _read_value_at_s3_key(key_path: str) -> str:
    """
    Read value at an s3 key

    Args:
        key_path (str): The S3 key path to read the value from.

    Returns:
        str: The value read from the specified S3 key path.
    """
    s3 = buckets.portal_resource

    wmg_config = CensusCubeConfig()
    wmg_config.load()

    s3obj = s3.Object(wmg_config.bucket, key_path)
    return s3obj.get()["Body"].read().decode("utf-8").strip()


def _should_reload_snapshot(
    *,
    snapshot_schema_version: str,
    explicit_snapshot_id_to_load: Optional[str] = None,
    snapshot_fs_root_path: Optional[str] = None,
) -> tuple[bool, str]:
    """
    Determine whether the in-memory snapshot should be reloaded and provide the id of the snapshot.

    Args:
        snapshot_schema_version (str): The version of the snapshot schema.
        explicit_snapshot_id_to_load (Optional[str]): The explicit snapshot id to load. Defaults to None.
        snapshot_fs_root_path (Optional[str]): The root path of the snapshot in the filesystem. Defaults to None.

    Returns:
        tuple[bool, str]: A pair of values. The first is a boolean indicating whether the in-memory snapshot should be reloaded. The second is the id of the snapshot that the in-memory data structure represents.
    """

    snapshot_id = explicit_snapshot_id_to_load or _get_latest_snapshot_id(
        snapshot_schema_version, snapshot_fs_root_path
    )

    if cached_snapshot is None:
        logger.info(f"Loading snapshot id: {snapshot_id}")
        return (True, snapshot_id)
    ######################### AN IMPORTANT NOTE #################################
    # 1. As of this writing on 01/12/2024, when the app is configured to read
    # the WMG snapshot from the local filesystem, the latest snapshot id will always
    # equal the snapshot id of the `cached_snapshot` object.
    #
    # This is because the scheme of downloading all the snapshots to the local filesystem
    # is done ONLY ONCE ON APP CONTAINER INITIALIZATION. That is, the code in the below
    # `elif` block will only execute if the app is configred to read from S3. If the app
    # is changed such that the local filesytem cache is updated by a another thread, then
    # the filesystem cache will behave like S3 in that updates to 'latest_snapshot_identifier'
    # file will eventually be reflected in the local disk while the app is still running. But
    # as of this writing, the WMG snapshot on local filesystem remains static throughout the
    # lifetime of a running application server process.
    #
    # 2. In the case of the application being configured to read from S3, note well that the
    # below `elif` condition will be come True if the value in 'latest_snapshot_identifier' file
    # is updated to be different from 'cached_snapshot.snapshot_identifer'. That is, it is possible
    # 'latest_snapshot_identifier' file could be modified to contain a snapshot_id that is older
    # than what is in 'cached_snapshot.snapshot_identifier'.
    #
    # This might be useful if we want to remove a corrupt wmg snapshot by simply
    # updating the latest_snapshot_identifer file to be an older snapshot_id without
    # requiring updating the explicit_snapshot_id_to_load in code and deploying.
    # Such an approach to rollback (updating the latest_snapshot_identifier on S3/filesystem)
    # should only be used in EXTREME emergencies. The normal rollback of updating the code's
    # config file and redeploying should be used in all other scenarios as this form of
    # rollback provides an audit log via git commit history.

    elif snapshot_id != cached_snapshot.snapshot_identifier:
        logger.info(
            f"Reloading snapshot. Detected snapshot id update from cached "
            f"snapshot id: {cached_snapshot.snapshot_identifier} to {snapshot_id}"
        )
        return (True, snapshot_id)
    else:
        logger.debug(f"No need to load snapshot. Using cached snapshot id: {cached_snapshot.snapshot_identifier}")
        return (False, snapshot_id)
