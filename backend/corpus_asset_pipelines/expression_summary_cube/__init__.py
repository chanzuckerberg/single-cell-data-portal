from backend.corpus_asset_pipelines.expression_summary_cube.job import create_expression_summary_cube


def run(corpus_path):
    create_expression_summary_cube(corpus_path)