"""
Drops and recreates all tables according to corpora_orm.py
"""

from sqlalchemy import create_engine

from backend.common.corpora_config import CorporaDbConfig
from backend.layers.persistence.orm import metadata


def create_db():
    engine = create_engine(CorporaDbConfig().database_uri)
    print("Dropping tables")
    metadata.drop_all(engine)
    print("Recreating tables")
    metadata.create_all(engine)


if __name__ == "__main__":
    create_db()
