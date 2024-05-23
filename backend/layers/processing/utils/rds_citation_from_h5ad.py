import re

from backend.common.utils.regex import DATASET_ID_REGEX


def rds_citation_from_h5ad(citation: str) -> str:
    return re.sub(f"{DATASET_ID_REGEX}\.h5ad", r"\g<dataset_id>.rds", citation)
