import logging
import time

from backend.cellguide.pipeline.source_collections import run as run_source_collections_pipeline

logging.basicConfig(level=logging.INFO)


def run_pipeline():
    output_directory = f"source_collections__{int(time.time())}"

    # Generate source data for each cell type
    run_source_collections_pipeline(output_directory=output_directory)


if __name__ == "__main__":
    run_pipeline()
