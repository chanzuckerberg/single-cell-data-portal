"""
Drops and recreates all tables according to corpora_orm.py
"""

import os
import sys

from sqlalchemy import create_engine

from backend.corpora.common.corpora_config import CorporaDbConfig
from backend.corpora.common.corpora_orm import Base


def create_db():
    engine = create_engine(CorporaDbConfig().database_uri)
    print("Dropping tables")
    Base.metadata.drop_all(engine)
    print("Recreating tables")
    Base.metadata.create_all(engine)


if __name__ == "__main__":
    create_db()
