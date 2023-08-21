import logging
import time

from backend.cellguide.pipeline.metadata import run as run_metadata_pipeline
from backend.cellguide.pipeline.ontology_tree import run as run_ontology_tree_pipeline

logging.basicConfig(level=logging.INFO)


def run_pipeline():
    output_directory = f"cellguide_pipeline_output__{int(time.time())}"

    # Run ontology tree pipeline
    ontology_tree = run_ontology_tree_pipeline(output_directory, get_tree_builder_only=True)

    # Generate cell guide cards, synonyms, and descriptions
    run_metadata_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)


if __name__ == "__main__":
    run_pipeline()
