import time

from backend.cellguide.pipeline.ontology_tree import run as run_ontology_tree_pipeline


def run_cellguide_pipeline():
    output_directory = f"cellguide_pipeline_output__{int(time.time())}"

    # Run ontology tree pipeline
    run_ontology_tree_pipeline(output_directory)
