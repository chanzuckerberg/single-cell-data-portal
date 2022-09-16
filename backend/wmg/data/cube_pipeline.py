import json
import logging
import pathlib
import sys
import time

from backend.corpora.common.utils.slack import (
    notify_slack,
    format_failed_batch_issue_slack_alert,
    gen_wmg_pipeline_success_message,
    gen_wmg_pipeline_failure_message,
)
from backend.corpus_asset_pipelines import integrated_corpus
from backend.corpus_asset_pipelines import summary_cubes

from backend.wmg.data.load_cube import upload_artifacts_to_s3, make_snapshot_active
from backend.wmg.data.transform import (
    generate_primary_filter_dimensions,
    get_cell_types_by_tissue,
    generate_cell_ordering,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def load_data_and_create_cube(
    path_to_h5ad_datasets: str,
    corpus_name: str = "corpus_group",
    snapshot_path=None,
    extract_data=True,
    validate_cube=True,
) -> (int, dict):
    """
    Function to copy H5AD datasets (from a preconfiugred s3 bucket) to the path given then,
    open, transform, normalize and concatenate them together as a tiledb object with a global gene index
    under the given corpus name.
    A cube of expression summary statistics (known as the expression summary cube) across genes (but queryable by
    the given dimensions/data) is then generated and stored under big-cube.
    A per-tissue mapping of cell ontologies is generated and the files are copied to s3 under a shared timestamp,.
    On success the least recent set of files are removed from s3
    """
    if not snapshot_path:
        snapshot_id = int(time.time())
        snapshot_path = f"{pathlib.Path().resolve()}/{snapshot_id}"
    corpus_path = f"{snapshot_path}/{corpus_name}"

    dataset_count = integrated_corpus.run(path_to_h5ad_datasets, corpus_path, extract_data)
    stats = summary_cubes.run(corpus_path, validate_cube)
    stats["dataset_count"] = dataset_count
    cell_type_by_tissue = get_cell_types_by_tissue(corpus_path)
    generate_cell_ordering(snapshot_path, cell_type_by_tissue)
    generate_primary_filter_dimensions(snapshot_path, corpus_name, snapshot_id)
    upload_artifacts_to_s3(snapshot_path, snapshot_id)
    if validate_cube:
        make_snapshot_active(snapshot_id)
        logger.info(f"Updated latest_snapshot_identifier in s3. Current snapshot id is {snapshot_id}")
        return snapshot_id, stats


def main():
    """
    To trigger the batch process in prod via the cli run
    AWS_PROFILE=single-cell-prod aws batch submit-job --job-name $JOB_NAME --job-queue dp-prod --job-definition dp-prod-prodstack-wmg-processing # noqa E501
    """
    # todo pass in validate_cubes as env arg
    try:
        snapshot_id, stats = load_data_and_create_cube("datasets", ".")
        pipeline_success_message = gen_wmg_pipeline_success_message(snapshot_id, stats)
        data = json.dumps(pipeline_success_message, indent=2)
        notify_slack(data)
    except Exception as e:
        logger.exception("Pipeline failed")
        pipeline_failure_message = gen_wmg_pipeline_failure_message(
            f"Issue with cube creation pipeline: {e}. See logs for more detail"
        )
        data = format_failed_batch_issue_slack_alert(pipeline_failure_message)
        notify_slack(data)


if __name__ == "__main__":
    main()
    sys.exit()
