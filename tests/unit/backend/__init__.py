# TODO(mbarrien): This should not be calling TestDatabase unconditionally upon import.

from tests.unit.fixtures.test_db import TestDatabase


testdb = TestDatabase()
testdb.create_db()
testdb.populate_test_data()
