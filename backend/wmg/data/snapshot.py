import json
import logging
import os
from dataclasses import dataclass
from typing import Dict, Optional

import pandas as pd
import tiledb
from pandas import DataFrame
from tiledb import Array

from backend.common.utils.s3_buckets import buckets
from backend.wmg.config import WmgConfig
from backend.wmg.data.schemas.corpus_schema import (
    DATASET_TO_GENE_IDS_NAME,
    FILTER_RELATIONSHIPS_NAME,
)
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import get_datasets_from_curation_api

# Snapshot data artifact file/dir names
CELL_TYPE_ORDERINGS_FILENAME = "cell_type_orderings.json"
PRIMARY_FILTER_DIMENSIONS_FILENAME = "primary_filter_dimensions.json"
EXPRESSION_SUMMARY_CUBE_NAME = "expression_summary"
EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME = "expression_summary_default"
EXPRESSION_SUMMARY_FMG_CUBE_NAME = "expression_summary_fmg"
CELL_COUNTS_CUBE_NAME = "cell_counts"
MARKER_GENES_CUBE_NAME = "marker_genes"
DATASET_TO_GENE_IDS_FILENAME = f"{DATASET_TO_GENE_IDS_NAME}.json"
FILTER_RELATIONSHIPS_FILENAME = f"{FILTER_RELATIONSHIPS_NAME}.json"

logger = logging.getLogger("wmg")


@dataclass
class WmgSnapshot:
    """
    All of the data artifacts the WMG API depends upon to perform its functions, versioned by "snapshot_identifier".
    These are read from data artifacts, per the relative file names, above.
    """

    snapshot_identifier: str

    # TileDB array containing expression summary statistics (expressed gene count, non-expressed mean,
    # etc.) aggregated by multiple cell metadata dimensions and genes. See the full schema at
    # backend/wmg/data/schemas/cube_schema.py.
    expression_summary_cube: Array

    # TileDB array containing expression summary statistics optimized for marker gene computation.
    # See the full schema at backend/wmg/data/schemas/expression_summary_fmg_cube_schema.py.
    expression_summary_fmg_cube: Array

    # TileDB array containing expression summary statistics optimized for querying with no
    # secondary filters selected.
    # See the full schema at backend/wmg/data/schemas/cube_schema_default.py.
    expression_summary_default_cube: Array

    # TileDB array containing the precomputed marker genes.
    # See the full schema at backend/wmg/data/schemas/marker_gene_cube_schema.py.
    marker_genes_cube: Array

    # TileDB array containing the total cell counts (expressed gene count, non-expressed mean, etc.) aggregated by
    # multiple cell metadata dimensions (but no gene dimension). See the full schema at
    # backend/wmg/data/schemas/cube_schema.py.
    cell_counts_cube: Array

    # Pandas DataFrame containing per-tissue ordering of cell types.
    # Columns are "tissue_ontology_term_id", "cell_type_ontology_term_id", "order"
    cell_type_orderings: DataFrame

    # precomputed list of ids for all gene and tissue ontology term ids per organism
    primary_filter_dimensions: Dict

    # dictionary of gene IDs mapped to dataset IDs
    dataset_to_gene_ids: Dict

    # precomputed filter relationships graph
    filter_relationships: Dict

    def __hash__(self):
        return hash(None)  # hash is not used for WmgSnapshot

    def build_dataset_metadata_dict(self):
        datasets = get_datasets_from_curation_api()
        dataset_dict = {}
        for dataset in datasets:
            dataset_id = dataset["dataset_id"]
            dataset_dict[dataset_id] = dict(
                id=dataset_id,
                label=dataset["title"],
                collection_id=dataset["collection_id"],
                collection_label=dataset["collection_name"],
            )
        self.dataset_dict = dataset_dict


# Cached data
cached_snapshot: Optional[WmgSnapshot] = None


def load_snapshot() -> WmgSnapshot:
    """
    Loads and caches the WMG snapshot. Reloads the snapshot data if the latest_snapshot_identifier S3 object has
    been updated.
    @return: WmgSnapshot object
    """
    global cached_snapshot
    new_snapshot_identifier = _update_latest_snapshot_identifier()
    if cached_snapshot is None and new_snapshot_identifier is not None:
        if os.environ.get("DEPLOYMENT_STAGE") == "test":
            cached_snapshot = _load_snapshot(new_snapshot_identifier)
        else:
            cached_snapshot = _download_and_load_snapshot(new_snapshot_identifier)
        cached_snapshot.build_dataset_metadata_dict()
    return cached_snapshot


def _load_snapshot(new_snapshot_identifier) -> WmgSnapshot:
    cell_type_orderings = _load_cell_type_order(new_snapshot_identifier)
    primary_filter_dimensions = _load_primary_filter_data(new_snapshot_identifier)
    dataset_to_gene_ids = _load_dataset_to_gene_ids_data(new_snapshot_identifier)
    filter_relationships = _load_filter_graph_data(new_snapshot_identifier)

    snapshot_base_uri = _build_snapshot_base_uri(new_snapshot_identifier)
    logger.info(f"Loading WMG snapshot at {snapshot_base_uri}")

    # TODO: Okay to keep TileDB arrays open indefinitely? Is it faster than re-opening each request?
    #  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell
    #  -data-portal/2134
    return WmgSnapshot(
        snapshot_identifier=new_snapshot_identifier,
        expression_summary_cube=_open_cube(f"{snapshot_base_uri}/{EXPRESSION_SUMMARY_CUBE_NAME}"),
        expression_summary_default_cube=_open_cube(f"{snapshot_base_uri}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}"),
        expression_summary_fmg_cube=_open_cube(f"{snapshot_base_uri}/{EXPRESSION_SUMMARY_FMG_CUBE_NAME}"),
        marker_genes_cube=_open_cube(f"{snapshot_base_uri}/{MARKER_GENES_CUBE_NAME}"),
        cell_counts_cube=_open_cube(f"{snapshot_base_uri}/{CELL_COUNTS_CUBE_NAME}"),
        cell_type_orderings=cell_type_orderings,
        primary_filter_dimensions=primary_filter_dimensions,
        dataset_to_gene_ids=dataset_to_gene_ids,
        filter_relationships=filter_relationships,
    )


def _download_and_load_snapshot(new_snapshot_identifier) -> WmgSnapshot:
    snapshot_base_uri = _build_snapshot_base_uri(new_snapshot_identifier)
    logger.info(f"Downloading WMG snapshot at {snapshot_base_uri}")

    local_snapshot_name = f"local_snapshot_{new_snapshot_identifier}"
    os.system(f"aws s3 sync {snapshot_base_uri} {local_snapshot_name}/")

    cell_type_orderings = _load_cell_type_order(local_snapshot_name)
    primary_filter_dimensions = _load_primary_filter_data(local_snapshot_name)
    dataset_to_gene_ids = _load_dataset_to_gene_ids_data(local_snapshot_name)
    filter_relationships = _load_filter_graph_data(local_snapshot_name)

    # TODO: Okay to keep TileDB arrays open indefinitely? Is it faster than re-opening each request?
    #  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell
    #  -data-portal/2134
    return WmgSnapshot(
        snapshot_identifier=new_snapshot_identifier,
        expression_summary_cube=_open_cube(f"{local_snapshot_name}/{EXPRESSION_SUMMARY_CUBE_NAME}"),
        expression_summary_default_cube=_open_cube(f"{local_snapshot_name}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}"),
        expression_summary_fmg_cube=_open_cube(f"{local_snapshot_name}/{EXPRESSION_SUMMARY_FMG_CUBE_NAME}"),
        marker_genes_cube=_open_cube(f"{local_snapshot_name}/{MARKER_GENES_CUBE_NAME}"),
        cell_counts_cube=_open_cube(f"{local_snapshot_name}/{CELL_COUNTS_CUBE_NAME}"),
        cell_type_orderings=cell_type_orderings,
        primary_filter_dimensions=primary_filter_dimensions,
        dataset_to_gene_ids=dataset_to_gene_ids,
        filter_relationships=filter_relationships,
    )


def _open_cube(cube_uri) -> Array:
    return tiledb.open(cube_uri, ctx=create_ctx(json.loads(WmgConfig().tiledb_config_overrides)))


def _load_cell_type_order(local_snapshot_name: str) -> DataFrame:
    return pd.read_json(f"{local_snapshot_name}/{CELL_TYPE_ORDERINGS_FILENAME}")


def _load_primary_filter_data(local_snapshot_name: str) -> Dict:
    with open(f"{local_snapshot_name}/{PRIMARY_FILTER_DIMENSIONS_FILENAME}", "r") as f:
        return json.load(f)


def _load_dataset_to_gene_ids_data(local_snapshot_name: str) -> Dict:
    with open(f"{local_snapshot_name}/{DATASET_TO_GENE_IDS_FILENAME}", "r") as f:
        return json.load(f)


def _load_filter_graph_data(local_snapshot_name: str) -> str:
    with open(f"{local_snapshot_name}/{FILTER_RELATIONSHIPS_FILENAME}", "r") as f:
        return json.load(f)


def _read_s3obj(relative_path: str) -> str:
    s3 = buckets.portal_resource
    wmg_config = WmgConfig()
    wmg_config.load()
    prefixed_relative_path = os.path.join(_build_data_path_prefix(), relative_path or "")
    s3obj = s3.Object(WmgConfig().bucket, prefixed_relative_path)
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


def _build_data_path_prefix():
    wmg_config = WmgConfig()
    rdev_prefix = os.environ.get("REMOTE_DEV_PREFIX")
    if rdev_prefix:
        return rdev_prefix.strip("/")
    elif "data_path_prefix" in wmg_config.config:
        return wmg_config.data_path_prefix
    else:
        return ""


def _build_snapshot_base_uri(snapshot_identifier: str):
    wmg_config = WmgConfig()
    data_path_prefix = _build_data_path_prefix()
    return os.path.join(
        "s3://",
        wmg_config.bucket,
        data_path_prefix,
        snapshot_identifier,
    )
