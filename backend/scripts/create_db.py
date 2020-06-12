import os
import sys

from sqlalchemy import create_engine

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from corpora.common.corpora_config import CorporaDbConfig
from corpora.common.corpora_orm import Base

conn_uri = CorporaDbConfig().database_uri
connection = f"{conn_uri[:conn_uri.index('@')]}@localhost:5432"  # TODO: clean up
engine = create_engine(connection)

print("Dropping tables")
Base.metadata.drop_all(engine)
print("Recreating tables")
Base.metadata.create_all(engine)
