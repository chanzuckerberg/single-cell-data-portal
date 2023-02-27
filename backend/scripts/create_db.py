"""
Drops and recreates all tables according to orm.py
"""


def current_db():
    from backend.layers.persistence.persistence import DatabaseProvider

    db_provider = DatabaseProvider()
    print("current db: Dropping tables")
    db_provider._drop_schema()
    print("current db: Recreating tables")
    db_provider._create_schema()


def create_db():
    current_db()


if __name__ == "__main__":
    create_db()
