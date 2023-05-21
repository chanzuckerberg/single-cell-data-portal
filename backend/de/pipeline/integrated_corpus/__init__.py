import tiledb

from backend.wmg.data.schemas.corpus_schema import create_tdb_integrated_corpus
from backend.wmg.pipeline.integrated_corpus.job import build_integrated_corpus, extract_datasets


def run(dataset_directory: list, corpus_path: str, extract_data: bool):
    """
    Copy relevant datasets and integrate cells into data corpus matrix
    """
    if not tiledb.VFS().is_dir(corpus_path):
        create_tdb_integrated_corpus(corpus_path)
    if extract_data:
        extract_datasets(dataset_directory)
    build_integrated_corpus(dataset_directory, corpus_path)
