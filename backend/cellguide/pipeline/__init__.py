import logging
import os
import shutil
import sys
import time
from glob import glob

from backend.cellguide.pipeline.canonical_marker_genes import run as run_canonical_marker_gene_pipeline
from backend.cellguide.pipeline.computational_marker_genes import run as run_computational_marker_gene_pipeline
from backend.cellguide.pipeline.config import CellGuideConfig
from backend.cellguide.pipeline.constants import GPT_OUTPUT_DIRECTORY_FOLDERNAME, GPT_SEO_OUTPUT_DIRECTORY_FOLDERNAME
from backend.cellguide.pipeline.gpt_descriptions import run as run_gpt_description_pipeline
from backend.cellguide.pipeline.metadata import run as run_metadata_pipeline
from backend.cellguide.pipeline.ontology_tree import run as run_ontology_tree_pipeline
from backend.cellguide.pipeline.providers.s3_provider import S3Provider
from backend.cellguide.pipeline.source_collections import run as run_source_collections_pipeline
from backend.cellguide.pipeline.utils import get_bucket_path, get_object_key
from backend.common.utils.cloudfront import create_invalidation_for_cellguide_data

logger = logging.getLogger(__name__)


def run_cellguide_pipeline():
    output_directory = str(int(time.time()))

    # delete any existing pipeline outputs
    cleanup(output_directory=output_directory)

    # Run ontology tree pipeline
    run_ontology_tree_pipeline(output_directory=output_directory)

    # Generate cell guide cards, synonyms, and descriptions
    run_metadata_pipeline(output_directory=output_directory)

    # Generate canonical marker genes from ASCT-B (HUBMAP)
    run_canonical_marker_gene_pipeline(output_directory=output_directory)

    # Generate computational marker genes from the CZI corpus
    run_computational_marker_gene_pipeline(output_directory=output_directory)

    # Generate source data for each cell type
    run_source_collections_pipeline(output_directory=output_directory)

    # Generate ChatGPT descriptions for any new cell ids
    run_gpt_description_pipeline(
        gpt_output_directory=GPT_OUTPUT_DIRECTORY_FOLDERNAME,
        gpt_seo_output_directory=GPT_SEO_OUTPUT_DIRECTORY_FOLDERNAME,
    )

    upload_cellguide_pipeline_output_to_s3(output_directory=output_directory)
    upload_gpt_descriptions_to_s3(
        gpt_output_directory=GPT_OUTPUT_DIRECTORY_FOLDERNAME,
        gpt_seo_output_directory=GPT_SEO_OUTPUT_DIRECTORY_FOLDERNAME,
    )

    # invalidate cloudfront distribution to reset cache
    create_invalidation_for_cellguide_data()

    # cleanup
    cleanup(output_directory=output_directory)


def upload_cellguide_pipeline_output_to_s3(*, output_directory: str):
    """
    If the pipeline is running in a deployed environment, then this function uploads
    the CellGuide snapshot to the corresponding environment's CellGuide data bucket.
    """

    s3_provider = S3Provider()
    bucket = CellGuideConfig().bucket

    bucket_path = get_bucket_path()
    if bucket_path is None:
        return

    logger.info(f"Uploading the pipeline output at {output_directory} to {bucket_path}")
    s3_provider.sync_directory(src_dir=output_directory, s3_uri=f"{bucket_path}{output_directory}")

    with open("latest_snapshot_identifier", "w") as file:
        file.write(output_directory)

    object_key = get_object_key(object="latest_snapshot_identifier")
    s3_provider.upload_file("latest_snapshot_identifier", bucket, object_key, {})

    # Create an empty file named 404
    with open("404", "w") as file:
        pass
    # Upload the empty 404 file to the bucket
    # this is used for custom cloudfront error handling
    s3_provider.upload_file("404", bucket, "404", {})


def upload_gpt_descriptions_to_s3(*, gpt_output_directory: str, gpt_seo_output_directory: str) -> None:
    bucket_path = get_bucket_path()

    s3_provider = S3Provider()
    for src_directory, dst_directory in zip(
        [gpt_output_directory, gpt_seo_output_directory],
        [GPT_OUTPUT_DIRECTORY_FOLDERNAME, GPT_SEO_OUTPUT_DIRECTORY_FOLDERNAME],
    ):
        # upload to s3
        s3_provider.sync_directory(src_dir=src_directory, s3_uri=f"{bucket_path}{dst_directory}/")

        num_descriptions = len(glob(f"{src_directory}/*.json"))
        logger.info(f"Uploaded {num_descriptions} GPT descriptions to {bucket_path}{dst_directory}/")


def cleanup(*, output_directory: str):
    logger.info(f"Cleaning up {output_directory} and other CellGuide pipeline outputs")
    if os.path.exists(output_directory):
        shutil.rmtree(output_directory)
    if os.path.exists(GPT_OUTPUT_DIRECTORY_FOLDERNAME):
        shutil.rmtree(GPT_OUTPUT_DIRECTORY_FOLDERNAME)
    if os.path.exists(GPT_SEO_OUTPUT_DIRECTORY_FOLDERNAME):
        shutil.rmtree(GPT_SEO_OUTPUT_DIRECTORY_FOLDERNAME)
    if os.path.exists("404"):
        os.remove("404")
    if os.path.exists("latest_snapshot_identifier"):
        os.remove("latest_snapshot_identifier")


if __name__ == "__main__":
    run_cellguide_pipeline()
    sys.exit()
