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
