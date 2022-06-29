from backend.corpus_asset_pipelines.integrated_corpus.job import build_integrated_corpus, extract_datasets


def run(dataset_directory: list, corpus_path: str, extract_data: bool):
    if extract_data:
        extract_datasets(dataset_directory)
    build_integrated_corpus(dataset_directory, corpus_path)
