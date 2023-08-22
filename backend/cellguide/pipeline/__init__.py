import logging
import os
import subprocess
import sys
import time

from backend.cellguide.pipeline.canonical_marker_genes import run as run_canonical_marker_gene_pipeline
from backend.cellguide.pipeline.computational_marker_genes import run as run_computational_marker_gene_pipeline
from backend.cellguide.pipeline.constants import CELL_GUIDE_DATA_BUCKET_PATH_PREFIX
from backend.cellguide.pipeline.metadata import run as run_metadata_pipeline
from backend.cellguide.pipeline.ontology_tree import run as run_ontology_tree_pipeline
from backend.cellguide.pipeline.source_collections import run as run_source_collections_pipeline

logger = logging.getLogger(__name__)


def run_cellguide_pipeline():
    output_directory = f"cellguide_pipeline_output__{int(time.time())}"

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


def upload_cellguide_pipeline_output_to_s3(*, output_directory: str):
    """
    If the pipeline is running in a deployed environment, then this function uploads
    the CellGuide snapshot to the corresponding environment's CellGuide data bucket.

    If the pipeline is running locally, this function uploads to the data bucket corresponding
    to the environment specified in CELLGUIDE_PIPELINE_TARGET_DEPLOYMENT
    """
    deployment_stage = os.getenv("DEPLOYMENT_STAGE")
    target_deployment = os.getenv("CELLGUIDE_PIPELINE_TARGET_DEPLOYMENT")

    if not deployment_stage and not target_deployment:
        logger.warning(
            f"Not uploading the pipeline output at {output_directory} to S3 as neither DEPLOYMENT_STAGE nor CELLGUIDE_PIPELINE_TARGET_DEPLOYMENT are set. Please set CELLGUIDE_PIPELINE_TARGET_DEPLOYMENT to one of dev, staging, or prod"
        )
        return

    if deployment_stage in ["dev", "staging", "prod"]:
        bucket_path = f"{CELL_GUIDE_DATA_BUCKET_PATH_PREFIX}{deployment_stage}/"
    elif target_deployment in ["dev", "staging", "prod"]:
        bucket_path = f"{CELL_GUIDE_DATA_BUCKET_PATH_PREFIX}{target_deployment}/"
    else:
        logger.warning(
            f"Invalid CELLGUIDE_PIPELINE_TARGET_DEPLOYMENT value: {target_deployment}. Please set CELLGUIDE_PIPELINE_TARGET_DEPLOYMENT to one of dev, staging, or prod"
        )
        return

    logger.info(f"Uploading the pipeline output at {output_directory} to {bucket_path}")
    sync_command = ["aws", "s3", "sync", f"{output_directory}/", bucket_path, "--quiet"]
    subprocess.run(sync_command)


if __name__ == "__main__":
    run_cellguide_pipeline()
    sys.exit()
