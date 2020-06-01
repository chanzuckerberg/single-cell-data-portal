import os
import sys

from sqlalchemy import create_engine

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from .corpora_config import CorporaDbConfig
from .corpora_orm import Base


class CorporaDatabase:
    def __init__(self):
        config = CorporaDbConfig()
        self.engine = create_engine(config.database_uri)

        print("Dropping tables")
        Base.metadata.drop_all(self.engine)
        print("Recreating tables")
        Base.metadata.create_all(self.engine)
        print([t[0] for t in self.engine.execute("SHOW TABLES;")])
