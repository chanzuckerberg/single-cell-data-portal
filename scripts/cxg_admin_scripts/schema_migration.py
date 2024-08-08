import json
from pathlib import Path

from backend.common.utils.json import CustomJSONEncoder
from backend.layers.processing.schema_migration import SchemaMigrate


def generate_report(ctx, execution_id: str, report_path: str, artifact_bucket: str):
    schema_migration = SchemaMigrate(ctx.obj["business_logic"], None)
    report = schema_migration.report(execution_id=execution_id, artifact_bucket=artifact_bucket, dry_run=True)
    report_file = Path(report_path).joinpath(f"{ctx.obj['deployment']}-{execution_id}.json")
    with open(report_file, "w") as f:
        json.dump(report, f, indent=4, sort_keys=True, cls=CustomJSONEncoder)
    print(f"Report saved to {report_file}")
