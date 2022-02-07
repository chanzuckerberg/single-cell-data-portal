"""
Drops and recreates all tables according to corpora_orm.py
"""

from sqlalchemy import create_engine

from backend.api.data_portal.config.app_config import DbConfig
from backend.api.data_portal.common.corpora_orm import Base


def create_db():
    engine = create_engine(DbConfig().database_uri)
    print("Dropping tables")
    Base.metadata.drop_all(engine)
    print("Recreating tables")
    Base.metadata.create_all(engine)


if __name__ == "__main__":
    create_db()
