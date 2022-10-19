import logging
import pathlib
import sys
import time

from backend.common.utils.result_notification import (
    notify_slack,
    format_failed_batch_issue_slack_alert,
    gen_wmg_pipeline_success_message,
    gen_wmg_pipeline_failure_message,
)
from backend.corpus_asset_pipelines import summary_cubes
from backend.corpus_asset_pipelines import integrated_corpus

from backend.wmg.data.load_cube import upload_artifacts_to_s3, make_snapshot_active
from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_CUBE_NAME, CELL_COUNTS_CUBE_NAME
from backend.wmg.data.transform import (
    generate_primary_filter_dimensions,
    get_cell_types_by_tissue,
    cell_type_ordering_create_file,
)
from backend.wmg.data.utils import (
    get_all_dataset_ids,
    get_expression_summary_cube_gene_count,
    get_cell_count_cube_count,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def load_data_and_create_cube(
    path_to_h5ad_datasets: str,
    corpus_name: str = "corpus_group",
    path=None,
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
    snapshot_id = int(time.time())
    snapshot_path = path if path else pathlib.Path().resolve()
    snapshot_path = f"{snapshot_path}/{snapshot_id}"
    corpus_path = f"{snapshot_path}/{corpus_name}"

    integrated_corpus.run(path_to_h5ad_datasets, corpus_path, extract_data)
    summary_cubes.run(corpus_path, validate_cube)
    stats = dict(
        dataset_count=len(get_all_dataset_ids(corpus_path)),
        gene_count=get_expression_summary_cube_gene_count(f"{corpus_path}/{EXPRESSION_SUMMARY_CUBE_NAME}"),
        cell_count=get_cell_count_cube_count(f"{corpus_path}/{CELL_COUNTS_CUBE_NAME}"),
    )
    cell_type_by_tissue = get_cell_types_by_tissue(corpus_path)
    cell_type_ordering_create_file(snapshot_path, cell_type_by_tissue)
    generate_primary_filter_dimensions(snapshot_path, corpus_name, snapshot_id)
    upload_artifacts_to_s3(snapshot_path, snapshot_id)
    if validate_cube:
        make_snapshot_active(snapshot_id)
        logger.info(f"Updated latest_snapshot_identifier in s3. Current snapshot id is {snapshot_id}")
    return snapshot_path, stats


def main():
    """
    To trigger the batch process in prod via the cli run
    AWS_PROFILE=single-cell-prod aws batch submit-job --job-name $JOB_NAME --job-queue dp-prod --job-definition dp-prod-prodstack-wmg-processing # noqa E501
    """
    # todo pass in validate_cubes as env arg
    try:
        snapshot_path, stats = load_data_and_create_cube("datasets", ".")
        success_message = gen_wmg_pipeline_success_message(snapshot_path, stats)
        notify_slack(success_message)
    except Exception as e:
        logger.exception("Pipeline failed")
        failure_message = format_failed_batch_issue_slack_alert(
            gen_wmg_pipeline_failure_message(
                f"Issue with WMG snapshot generation pipeline: {e}. See logs for more detail."
            )
        )
        notify_slack(failure_message)


if __name__ == "__main__":
    main()
    sys.exit()
