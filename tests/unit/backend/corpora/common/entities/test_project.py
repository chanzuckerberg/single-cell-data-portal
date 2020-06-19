import os
import sys
import unittest

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..", ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.entities.project import Project


class TestProject(unittest.TestCase):
    def setUp(self):
        self.id = "f753497f-f8a5-48b5-b344-de9ad5a4354d"
        self.status = "LIVE"

    @unittest.skipIf(True, "Only runnable on local dev env via bastion tunnel. Comment this line to run.")
    def test_get(self):
        # TODO: write tests
        project = Project.get((self.id, self.status))
        print(project.validation_state)
        print(project.dataset_ids)
        print(project.contributors[0]['name'])
