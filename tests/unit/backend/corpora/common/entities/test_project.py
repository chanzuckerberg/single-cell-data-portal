import os
import sys
import unittest

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", "..", ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.entities.project import Project
from backend.corpora.common.utils.exceptions import CorporaException


class TestProject(unittest.TestCase):
    def setUp(self):
        self.uuid = "f753497f-f8a5-48b5-b344-de9ad5a4354d"
        self.status = "LIVE"

    @unittest.skipIf(True, "Only runnable on local dev env via bastion tunnel. Comment this line to run.")
    def test_get__ok(self):
        key = (self.uuid, self.status)

        project = Project.get(key)

        self.assertEqual(project.name, "Single-cell gene expression profiling of SARS-CoV-2 infected human cell lines")
        self.assertEqual(project.contributors[0]['name'], "Wyler Emanuel")

    @unittest.skipIf(True, "Only runnable on local dev env via bastion tunnel. Comment this line to run.")
    def test_get__dne(self):
        key = ("bad", self.status)

        with self.assertRaises(CorporaException):
            Project.get(key)

    @unittest.skipIf(True, "Only runnable on local dev env via bastion tunnel. Comment this line to run.")
    def test_get__invalid_status(self):
        key = (self.uuid, "bad")

        with self.assertRaises(CorporaException):
            Project.get(key)
