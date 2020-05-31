import unittest
from unittest import mock

from backend.corpora.common.utils.db_utils import DbUtils


class TestDbUtils(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.db = DbUtils()
        cls.db.create()

    def test_parse_multivalue(self):
        self.assertEqual(DbUtils._parse_multivalue(""), [])
        self.assertEqual(DbUtils._parse_multivalue(None), [])
        self.assertEqual(DbUtils._parse_multivalue("hello,world"), ["hello", "world"])

    @mock.patch("sqlalchemy.schema.MetaData.create_all")
    @mock.patch("sqlalchemy.schema.MetaData.drop_all")
    def test_create(self, mock_drop_all, mock_create_all):
        with self.subTest("Operation allowed against test db"):
            self.db.create()

            self.assertTrue(mock_drop_all.called)
            self.assertTrue(mock_create_all.called)

        with self.subTest("Throws EnvironmentError if not test db"):
            with mock.patch("backend.corpora.common.utils.db_utils.DbUtils._is_test_db") as mock_is_test_db:
                mock_is_test_db.return_value = False
                with self.assertRaises(EnvironmentError):
                    self.db.create()
