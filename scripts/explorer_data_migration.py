from urllib.parse import urlparse

import os
import sys

from backend.corpora.common.entities import Dataset
from backend.corpora.common.utils.db_session import db_session_manager, DBSessionMaker

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.corpora_config import CorporaDbConfig

def get_database_uri() -> str:
    uri = urlparse(CorporaDbConfig().database_uri)
    uri = uri._replace(netloc="@".join([uri[1].split("@")[0], "localhost:5432"]))
    return uri.geturl()

DBSessionMaker(get_database_uri())


def set_explorer_url_from_deployment_directory():
    with db_session_manager() as session:
        datasets = Dataset.list(session)
        for dataset in datasets:
            dataset.update(explorer_url=dataset.deployment_directories[0].url)

def set_deployment_directory_from_explorer_url():
    with db_session_manager() as session:
        datasets = Dataset.list(session)
        for dataset in datasets:
            dataset.update(deployment_directories=[dataset.explorer_url])

def set_explorer_s3_uri_from_explorer_url():
    pass


