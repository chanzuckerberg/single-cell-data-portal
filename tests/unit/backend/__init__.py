# TODO(mbarrien): This should not be calling TestDatabase unconditionally upon import.

import os
from ..backend.corpora.fixtures.database import TestDatabase


# TODO(jgadling): This env var workarond is a temporary fix - we're investigating a nicer refactor of unit tests.
if not os.getenv("SKIP_DB_RELOAD"):
    testdb = TestDatabase()
    testdb.create_db()
    testdb.populate_test_data()
