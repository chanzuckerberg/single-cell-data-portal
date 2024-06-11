import unittest

from backend.common.doi import clean_doi


class TestDoi(unittest.TestCase):
    def test__clean_doi(self):
        test_cases = [
            ("10.1016/j.cell.2019.11.025", "10.1016/j.cell.2019.11.025"),
            ("DOI: 10.1016/j.cell.2019.11.025.", "10.1016/j.cell.2019.11.025"),
            (" DOI: 10.1016/j.cell.2019.11.025 ", "10.1016/j.cell.2019.11.025"),
            ("10.1016/j.cell.2019.11.025. ", "10.1016/j.cell.2019.11.025"),
            ("", ""),
        ]

        for doi, expected in test_cases:
            with self.subTest(doi=doi):
                self.assertEqual(clean_doi(doi), expected)
