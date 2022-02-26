from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager


def get_s3_uris():
    with db_session_manager() as session:
        dataset = Dataset.get(session, include_tombstones=True)