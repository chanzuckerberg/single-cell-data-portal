import json
from pathlib import Path

from backend.common.utils.json import CustomJSONEncoder
from backend.layers.common.entities import CollectionVersionId, DatasetId
from backend.layers.processing.schema_migration import SchemaMigrate


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


def generate_report(ctx, execution_id: str, report_path: str, artifact_bucket: str):
    schema_migration = SchemaMigrate(ctx.obj["business_logic"], None)
    report = schema_migration.report(execution_id=execution_id, artifact_bucket=artifact_bucket, dry_run=True)
    report_file = Path(report_path).joinpath(f"{ctx.obj['deployment']}-{execution_id}.json")
    with open(report_file, "w") as f:
        json.dump(report, f, indent=4, sort_keys=True, cls=CustomJSONEncoder)
