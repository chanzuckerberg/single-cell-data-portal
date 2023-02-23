import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa


def get_public_dataset_details():
    # id, name, organism, tissue, assay, sex, cell_count, explorer_url, and S3 uris
    pass
