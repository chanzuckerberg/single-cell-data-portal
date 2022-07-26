"""
Recreates database.
"""

from backend.corpora.common.entities.tiledb_data import TileDBData

def create_db():
    print("Recreating db")
    TileDBData.init_db()


if __name__ == "__main__":
    create_db()
