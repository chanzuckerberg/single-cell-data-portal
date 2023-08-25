import logging
import os
import subprocess
import sys
import time

from backend.cellguide.pipeline.canonical_marker_genes import run as run_canonical_marker_gene_pipeline
from backend.cellguide.pipeline.computational_marker_genes import run as run_computational_marker_gene_pipeline
from backend.cellguide.pipeline.config import CellGuideConfig
from backend.cellguide.pipeline.metadata import run as run_metadata_pipeline
from backend.cellguide.pipeline.ontology_tree import run as run_ontology_tree_pipeline
from backend.cellguide.pipeline.source_collections import run as run_source_collections_pipeline
from backend.common.utils.cloudfront import create_invalidation_for_cellguide_data

logger = logging.getLogger(__name__)


def run_cellguide_pipeline():
    output_directory = str(int(time.time()))

    # Run ontology tree pipeline
    run_ontology_tree_pipeline(output_directory=output_directory)

    # Generate cell guide cards, synonyms, and descriptions
    run_metadata_pipeline(output_directory=output_directory)

    # Generate canonical marker genes from ASCT-B (HUBMAP)
    run_canonical_marker_gene_pipeline(output_directory=output_directory)

    # Generate source data for each cell type
    run_source_collections_pipeline(output_directory=output_directory)

    # Generate computational marker genes from the CZI corpus
    run_computational_marker_gene_pipeline(output_directory=output_directory)

    upload_cellguide_pipeline_output_to_s3(output_directory=output_directory)

    # invalidate cloudfront distribution to reset cache
    create_invalidation_for_cellguide_data()


def upload_cellguide_pipeline_output_to_s3(*, output_directory: str):
    """
    If the pipeline is running in a deployed environment, then this function uploads
    the CellGuide snapshot to the corresponding environment's CellGuide data bucket.
    """
    deployment_stage = os.getenv("DEPLOYMENT_STAGE")

    if not deployment_stage:
        logger.warning(f"Not uploading the pipeline output at {output_directory} to S3 as DEPLOYMENT_STAGE is not set.")
        return

    if deployment_stage in ["dev", "staging", "prod"]:
        bucket = CellGuideConfig().bucket
        bucket_path = f"s3://{bucket}/"
    else:
        logger.warning(
            f"Invalid DEPLOYMENT_STAGE value: {deployment_stage}. Please set DEPLOYMENT_STAGE to one of dev, staging, or prod"
        )
        return

    logger.info(f"Uploading the pipeline output at {output_directory} to {bucket_path}")
    sync_command = ["aws", "s3", "sync", f"{output_directory}/", f"{bucket_path}{output_directory}", "--quiet"]

    with open("latest_snapshot_identifier", "w") as file:
        file.write(output_directory)

    copy_command = [
        "aws",
        "s3",
        "cp",
        "latest_snapshot_identifier",
        f"{bucket_path}latest_snapshot_identifier",
        "--quiet",
    ]
    subprocess.run(copy_command)

    subprocess.run(sync_command)

    # Create an empty file named 404
    with open("404", "w") as file:
        pass

    # Upload the empty 404 file to the bucket
    # this is used for custom cloudfront error handling
    upload_404_command = ["aws", "s3", "cp", "404", f"{bucket_path}404", "--quiet"]
    subprocess.run(upload_404_command)


if __name__ == "__main__":
    run_cellguide_pipeline()
    sys.exit()
