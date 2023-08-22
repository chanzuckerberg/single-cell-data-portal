import logging
import time

from backend.cellguide.pipeline.computational_marker_genes import run as run_computational_marker_gene_pipeline
from backend.cellguide.pipeline.ontology_tree import run as run_ontology_tree_pipeline

logging.basicConfig(level=logging.INFO)


def run_pipeline():
    output_directory = f"cellguide_pipeline_output__{int(time.time())}"

    # Run ontology tree pipeline
    ontology_tree = run_ontology_tree_pipeline(output_directory, get_tree_builder_only=True)

    # Generate computational marker genes from the CZI corpus
    run_computational_marker_gene_pipeline(output_directory=output_directory, ontology_tree=ontology_tree)


if __name__ == "__main__":
    run_pipeline()
