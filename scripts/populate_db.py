#!/usr/bin/env python

import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from tests.unit.fixtures.test_db import TestDatabase

if __name__ == "__main__":
    testdb = TestDatabase(real_data=True)
    testdb.create_db()
    testdb.populate_test_data()
