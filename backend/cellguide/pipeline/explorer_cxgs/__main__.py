import logging
import time

from backend.cellguide.pipeline.explorer_cxgs import run as run_valid_explorer_cxgs_pipeline

logging.basicConfig(level=logging.INFO)


def run_pipeline():
    output_directory = f"valid_explorer_cxgs__{int(time.time())}"

    # Generate cell guide cards, synonyms, and descriptions
    run_valid_explorer_cxgs_pipeline(output_directory=output_directory)


if __name__ == "__main__":
    run_pipeline()
