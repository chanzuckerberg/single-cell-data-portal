import logging
import time

from backend.cellguide.pipeline.canonical_marker_genes import run as run_canonical_marker_gene_pipeline

logging.basicConfig(level=logging.INFO)


def run_pipeline():
    output_directory = f"cellguide_pipeline_output__{int(time.time())}"

    # Generate computational marker genes from the CZI corpus
    run_canonical_marker_gene_pipeline(output_directory=output_directory)


if __name__ == "__main__":
    run_pipeline()
