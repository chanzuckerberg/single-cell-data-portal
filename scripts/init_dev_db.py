#!/usr/bin/env python

import os
import sys

from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, create_database

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.corpora_config import CorporaDbConfig

if __name__ == "__main__":
    engine = create_engine(CorporaDbConfig().database_uri)

    if not database_exists(engine.url):
        print("Database does not exist, creating database")
        create_database(engine.url)
        # TODO(mbarrien): HACK HACK HACK! TestDatabase drops database upon import.
        # Import inside function to avoid that undesired behavior.
        from tests.unit.backend.corpora.fixtures.database import TestDatabase
        TestDatabase(real_data=True)
