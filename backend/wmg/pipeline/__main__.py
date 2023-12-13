import logging

from backend.wmg.pipeline import run_pipeline

logging.basicConfig(level=logging.INFO)

if __name__ == "__main__":
    run_pipeline()
