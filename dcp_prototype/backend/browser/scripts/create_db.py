import os
import sys

from sqlalchemy import create_engine

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.code.config.db_config import BrowserDbConfig
from browser.code.common import browser_orm as orm


if __name__ == "__main__":
    engine = create_engine(BrowserDbConfig().database_uri)
    print("Dropping tables")
    orm.Base.metadata.drop_all(engine)
    print("Recreating tables")
    orm.Base.metadata.create_all(engine)
    print([t[0] for t in engine.execute("SHOW TABLES;")])
