import json
from pathlib import Path

from backend.layers.common.entities import CollectionVersionId, DatasetId


def rollback_dataset(ctx, report_path: Path):
    with report_path.open("r") as f:
        report = json.load(f)
    for entry in report:
        if entry["rollback"] is True:
            collection_version_id = entry["collection_version_id"]
            dataset_id = entry["dataset_id"]
            ctx.obj["business_logic"].restore_previous_dataset_version(
                CollectionVersionId(collection_version_id), DatasetId(dataset_id)
            )


def rollback_collections_to_schema_version(ctx, schema_version: str):
    collection_map = ctx.obj["business_logic"].get_collection_map_to_latest_published_version_by_schema(schema_version)
    ctx.obj["business_logic"].set_collection_versions_as_canonical(collection_map)
