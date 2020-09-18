#!/usr/bin/env python

import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from tests.unit.backend.corpora.fixtures.database import TestDatabase

if __name__ == "__main__":
    TestDatabase(real_data=True)
