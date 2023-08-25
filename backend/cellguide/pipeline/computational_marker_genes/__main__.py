import logging
import time

from backend.cellguide.pipeline.computational_marker_genes import run as run_computational_marker_gene_pipeline

logging.basicConfig(level=logging.INFO)


def run_pipeline():
    output_directory = f"computational_marker_genes__{int(time.time())}"

    # Generate computational marker genes from the CZI corpus
    run_computational_marker_gene_pipeline(output_directory=output_directory)


if __name__ == "__main__":
    run_pipeline()
