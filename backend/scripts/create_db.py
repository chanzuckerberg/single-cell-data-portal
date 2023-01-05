"""
Drops and recreates all tables according to orm.py
"""
from backend.layers.persistence.persistence import DatabaseProvider


def create_db():
    db_provider = DatabaseProvider()
    print("Dropping tables")
    db_provider._drop_schema()
    print("Recreating tables")
    db_provider._create_schema()


if __name__ == "__main__":
    create_db()
