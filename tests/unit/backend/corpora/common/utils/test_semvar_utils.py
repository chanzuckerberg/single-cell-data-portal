import unittest

from backend.corpora.common.utils.semvar_utils import validate_version_str


class TestSemvarUtils(unittest.TestCase):
    def test__release_only_properly_formatted_semvar__returns_true(self):
        self.assertTrue(validate_version_str("0.17.34"))

    def test__release_only_misformatted_version__returns_false(self):
        self.assertFalse(validate_version_str("0.17"))

    def test__prerelease_properly_formatted_semvar__returns_true(self):
        self.assertTrue(validate_version_str("0.17.34-rc.1", release_only=False))

    def test__release_only_prerelease_properly_formatted_semvar__returns_false(self):
        self.assertFalse(validate_version_str("0.17.34-rc.1", release_only=True))

    def test__prerelease_release_version__returns_true(self):
        self.assertTrue(validate_version_str("0.17.34", release_only=False))