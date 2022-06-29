import sys

from backend.corpus_asset_pipelines.expression_summary_cube.job import create_expression_summary_cube
from backend.wmg.data.validation.validation import Validation


def run(corpus_path, validate_cube):
    create_expression_summary_cube(corpus_path)
    if validate_cube:
        if Validation(corpus_path).validate_cube() is False:
            sys.exit("Exiting due to cube validation failure")
