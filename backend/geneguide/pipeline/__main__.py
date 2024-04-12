import logging

from backend.cellguide.pipeline import run_cellguide_pipeline

logging.basicConfig(level=logging.INFO)

if __name__ == "__main__":
    run_cellguide_pipeline()
