import os
import sys
import unittest

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..", ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

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

    def test__get__dne(self):
        key = ("bad", self.status)

        self.assertEqual(Project.get(key), None)

    def test__get__invalid_status(self):
        key = (self.uuid, "bad")

        with self.assertRaises(CorporaException) as context:
            Project.get(key)
        self.assertEqual("Invalid status bad. Status must be one of ['LIVE', 'EDIT'].", str(context.exception))
