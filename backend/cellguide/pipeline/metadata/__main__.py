import logging
import time

from backend.cellguide.pipeline.metadata import run as run_metadata_pipeline

logging.basicConfig(level=logging.INFO)


def run_pipeline():
    output_directory = f"metadata__{int(time.time())}"

    # Generate cell guide cards, synonyms, and descriptions
    run_metadata_pipeline(output_directory=output_directory)


if __name__ == "__main__":
    run_pipeline()
