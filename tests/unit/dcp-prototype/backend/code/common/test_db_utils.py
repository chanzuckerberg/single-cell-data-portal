import json
import unittest
from unittest import mock

from dcp_prototype.backend.code.common.db_utils import DbUtils
from dcp_prototype.backend.scripts.load_artifact import load_from_artifact
from dcp_prototype.backend.scripts.mock import mock_data


class TestDbUtils(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path_to_fixtures = "tests/unit/dcp-prototype/backend/fixtures"
        test_artifact_file = f"{path_to_fixtures}/test_artifact.json"
        with open(test_artifact_file, "r") as f:
            cls.test_artifact = json.loads(f.read())

        cls.db = DbUtils()
        cls.db.create()
        load_from_artifact(session=cls.db.session, path_to_file=test_artifact_file)

    def setUp(self):
        self.test_project_id = "cc95ff89-2e68-4a08-a234-480eca21ce79"
        self.test_file_id = "ba17c4674f7356782757e895e8768f61"
        self.bad_id = "DNE"

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
            with mock.patch("dcp_prototype.backend.code.common.db_utils.DbUtils._is_test_db") as mock_is_test_db:
                mock_is_test_db.return_value = False
                with self.assertRaises(EnvironmentError):
                    self.db.create()

    def test_query_projects(self):
        projects = self.db.query_projects()

        self.assertEqual(len(projects), 3)

        self.assertEqual(projects[0]["id"], "48c1a1bc-31bd-4420-9d28-9e4e6289eb16")
        self.assertEqual(projects[0]["title"], "Single-cell RNA-seq of mouse cerebral cortex")
        self.assertEqual(projects[0]["assays"], ["Fluidigm C1-based library preparation"])
        self.assertEqual(projects[0]["organs"], ["CA1 field of hippocampus", "somatosensory cortex"])
        self.assertEqual(projects[0]["species"], ["Mus musculus"])
        self.assertEqual(projects[0]["cell_count"], mock_data[projects[0]["id"]]["cell_count"])

    def test_query_project(self):
        with self.subTest("Project exists"):
            project = self.db.query_project(self.test_project_id)

            self.assertEqual(project["id"], "cc95ff89-2e68-4a08-a234-480eca21ce79")
            self.assertEqual(project["title"], "1M Immune Cells")
            self.assertEqual(project["label"], mock_data[project["id"]]["label"])
            self.assertEqual(project["assays"], ["10X v2 sequencing"])
            self.assertEqual(project["organs"], ["blood", "immune system"])
            self.assertEqual(project["species"], ["Homo sapiens"])
            self.assertEqual(project["contributors"][0]["name"], "Aviv Regev")
            self.assertEqual(project["contributors"][0]["institution"], "Broad Institute")
            self.assertEqual(len(project["contributors"]), 14)
            self.assertEqual(project["description"], mock_data[project["id"]]["description"])
            self.assertEqual(project["biosample_categories"], ["primary tissue"])
            self.assertEqual(project["development_stages"], ["human adult stage", "postpartum"])
            self.assertEqual(project["diseases"], [])
            self.assertEqual(project["cell_isolation_methods"], ["10X v2 sequencing"])
            self.assertEqual(
                project["cell_types"], ["cord blood hematopoietic stem cell", "bone marrow hematopoietic cell"]
            )
            self.assertEqual(project["cell_count"], mock_data[project["id"]]["cell_count"])
            self.assertEqual(project["paired_end"], ["False"])
            self.assertEqual(project["nucleic_acid_sources"], ["single cell"])
            self.assertEqual(project["input_nucleic_acid_molecules"], ["polyA RNA"])
            self.assertEqual(project["publication_title"], "")
            self.assertEqual(project["publication_doi"], "")
            self.assertEqual(project["cxg_enabled"], mock_data[project["id"]]["cxg_enabled"])

        with self.subTest("Project DNE"):
            result = self.db.query_project(self.bad_id)
            self.assertIsNone(result)

    def test_query_file(self):
        with self.subTest("File exists"):
            file = self.db.query_file(self.test_file_id)

            self.assertEqual(file["id"], self.test_file_id)
            self.assertEqual(file["project_id"], "cc95ff89-2e68-4a08-a234-480eca21ce79")
            self.assertEqual(file["filename"], "matrix.loom")
            self.assertEqual(file["file_format"], "LOOM")
            self.assertEqual(file["file_size"], 3215269107)
            self.assertEqual(file["file_type"], "EXPRESSION_MATRIX")
            self.assertEqual(file["s3_uri"], "s3://demo-matrices-scd/1M Immune Cells/matrix.loom")

        with self.subTest("File DNE"):
            result = self.db.query_file(self.bad_id)
            self.assertIsNone(result)

    def test_query_project_assays(self):
        assays = self.db.query_project_assays(self.test_project_id)
        self.assertEqual(assays, ["10X v2 sequencing"])

    def test_query_project_organs(self):
        organs = self.db.query_project_organs(self.test_project_id)
        self.assertEqual(organs, ["blood", "immune system"])

    def test_query_project_species(self):
        species = self.db.query_project_species(self.test_project_id)
        self.assertEqual(species, ["Homo sapiens"])

    def test_query_project_contributors(self):
        contributors = self.db.query_project_contributors(self.test_project_id)

        self.assertEqual(len(contributors), 14)
        self.assertEqual(contributors[0]["name"], "Aviv Regev")
        self.assertEqual(contributors[0]["institution"], "Broad Institute")

    def test_query_downloadable_project_files(self):
        with self.subTest("Files found"):
            files = self.db.query_downloadable_project_files(self.test_project_id)

            self.assertEqual(len(files), 1)
            self.assertEqual(files[0]["id"], "ba17c4674f7356782757e895e8768f61")
            self.assertEqual(files[0]["filename"], "matrix.loom")
            self.assertEqual(files[0]["file_format"], "LOOM")
            self.assertEqual(files[0]["file_type"], "EXPRESSION_MATRIX")
            self.assertEqual(files[0]["file_size"], 3215269107)

        with self.subTest("No files found"):
            files = self.db.query_downloadable_project_files(self.bad_id)
            self.assertEqual(files, [])
