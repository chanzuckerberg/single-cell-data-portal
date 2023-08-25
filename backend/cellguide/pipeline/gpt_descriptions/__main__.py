import logging
import time

from backend.cellguide.pipeline.gpt_descriptions import run as run_gpt_description_pipeline

logging.basicConfig(level=logging.INFO)


def run_pipeline():
    output_directory = f"gpt_descriptions__{int(time.time())}"

    # Generate computational marker genes from the CZI corpus
    run_gpt_description_pipeline(output_directory=output_directory)


if __name__ == "__main__":
    run_pipeline()
