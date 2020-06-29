"""
Drops and recreates all tables according to corpora_orm.py
"""

import os
import sys

from sqlalchemy import create_engine

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from corpora.common.corpora_config import CorporaDbConfig
from corpora.common.corpora_orm import Base


def create_db():
    engine = create_engine(CorporaDbConfig().database_uri)
    print("Dropping tables")
    Base.metadata.drop_all(engine)
    print("Recreating tables")
    Base.metadata.create_all(engine)


if __name__ == "__main__":
    create_db()
