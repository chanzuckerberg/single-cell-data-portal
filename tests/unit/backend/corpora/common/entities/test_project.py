import unittest
from backend.corpora.common.entities.project import Project
from backend.corpora.common.utils.exceptions import CorporaException
from backend.corpora.common.corpora_orm import ProjectStatus


class TestProject(unittest.TestCase):
    def setUp(self):
        self.uuid = "test_project_id"
        self.status = ProjectStatus.LIVE.name

    def test__get__ok(self):
        key = (self.uuid, self.status)

        project = Project.get(key)

        self.assertEqual(project.name, "test_project")

    def test__get__does_not_exist(self):
        non_existent_key = ("non_existent_id", self.status)

        self.assertEqual(Project.get(non_existent_key), None)

    def test__get__invalid_status(self):
        invalid_status_key = (self.uuid, "invalid_status")

        with self.assertRaises(CorporaException) as context:
            Project.get(invalid_status_key)

        self.assertEqual(
            "Invalid status invalid_status. Status must be one of ['LIVE', 'EDIT'].", str(context.exception)
        )
