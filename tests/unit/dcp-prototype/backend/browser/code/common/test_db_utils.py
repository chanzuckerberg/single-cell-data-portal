import unittest
from unittest import mock

from dcp_prototype.backend.browser.code.common.db_utils import DbUtils


class TestDbUtils(unittest.TestCase):
    def setUp(self):
        self.project_id = "test_id"
        self.file_id = "file_id"

    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.query_project_species")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.query_project_organs")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.query_project_assays")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils._query")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.__init__")
    def test_query_projects(
        self, mock_init, mock_query, mock_query_project_assays, mock_query_project_organs, mock_query_projet_species
    ):
        mock_init.return_value = None
        mock_query.return_value = [mock.Mock()]

        db = DbUtils()
        projects = db.query_projects()
        self.assertEqual(len(projects), 1)
        self.assertIn("id", projects[0])
        self.assertIn("title", projects[0])
        self.assertIn("assays", projects[0])
        self.assertIn("organs", projects[0])
        self.assertIn("species", projects[0])
        self.assertIn("cell_count", projects[0])

    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.query_project_contributors")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.query_project_species")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.query_project_organs")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.query_project_assays")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils._get")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.__init__")
    def test_query_project(
        self,
        mock_init,
        mock_get,
        mock_query_project_assays,
        mock_query_project_organs,
        mock_query_projet_species,
        mock_query_project_contributors,
    ):
        mock_init.return_value = None
        db = DbUtils()

        with self.subTest("Project exists"):
            mock_get.return_value = mock.Mock()
            result = db.query_project(self.project_id)
            self.assertIn("id", result)
            self.assertIn("cxg_enabled", result)

        with self.subTest("Project DNE"):
            mock_get.return_value = None
            result = db.query_project(self.project_id)
            self.assertIsNone(result)

    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils._get")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.__init__")
    def test_query_file(self, mock_init, mock_get):
        mock_init.return_value = None
        db = DbUtils()

        with self.subTest("File exists"):
            mock_get.return_value = mock.Mock()
            result = db.query_file(self.file_id)
            self.assertIsNotNone(result)

        with self.subTest("File DNE"):
            mock_get.return_value = None
            result = db.query_file(self.file_id)
            self.assertIsNone(result)

    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils._query")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.__init__")
    def test_query_project_assays(self, mock_init, mock_query):
        mock_init.return_value = None
        mock_query.return_value = [mock.Mock()]
        db = DbUtils()

        result = db.query_project_assays(self.project_id)
        self.assertEqual(len(mock_query.return_value), len(result))

    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils._query")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.__init__")
    def test_query_project_organs(self, mock_init, mock_query):
        mock_init.return_value = None
        mock_query.return_value = [mock.Mock()]
        db = DbUtils()

        result = db.query_project_organs(self.project_id)
        self.assertEqual(len(mock_query.return_value), len(result))

    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils._query")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.__init__")
    def test_query_project_species(self, mock_init, mock_query):
        mock_init.return_value = None
        mock_query.return_value = [mock.Mock()]
        db = DbUtils()

        result = db.query_project_species(self.project_id)
        self.assertEqual(len(mock_query.return_value), len(result))

    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils._query")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.__init__")
    def test_query_project_contributors(self, mock_init, mock_query):
        mock_init.return_value = None
        mock_query.return_value = [mock.Mock()]
        db = DbUtils()

        result = db.query_project_contributors(self.project_id)
        self.assertEqual(len(mock_query.return_value), len(result))

    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils._query")
    @mock.patch("dcp_prototype.backend.browser.code.common.db_utils.DbUtils.__init__")
    def test_query_downloadable_project_files(self, mock_init, mock_query):
        mock_init.return_value = None
        mock_query.return_value = [mock.Mock()]

        db = DbUtils()

        result = db.query_downloadable_project_files(self.project_id)
        self.assertEqual(len(mock_query.return_value), len(result))
        self.assertIn("id", result[0])
        self.assertIn("filename", result[0])
