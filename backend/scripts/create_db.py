"""
Drops and recreates all tables according to corpora_orm.py and orm.py
"""


def legacy_db():
    from sqlalchemy import create_engine

    from backend.common.corpora_config import CorporaDbConfig
    from backend.common.corpora_orm import Base

    engine = create_engine(CorporaDbConfig().database_uri)
    print("legacy db: Dropping tables")
    Base.metadata.drop_all(engine)
    print("legacy db: Recreating tables")
    Base.metadata.create_all(engine)


def current_db():
    from backend.layers.persistence.persistence import DatabaseProvider

    db_provider = DatabaseProvider()
    print("current db: Dropping tables")
    db_provider._drop_schema()
    print("current db: Recreating tables")
    db_provider._create_schema()


def create_db():
    legacy_db()
    current_db()


if __name__ == "__main__":
    create_db()
