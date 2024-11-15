import unittest

from backend.common.doi import clean_doi


class TestDoi(unittest.TestCase):
    def test__clean_doi(self):
        test_cases = [
            ("10.1016/j.cell.2019.11.025", "10.1016/j.cell.2019.11.025"),
            ("DOI: 10.1016/j.cell.2019.11.025.", "10.1016/j.cell.2019.11.025"),
            (" DOI: 10.1016/j.cell.2019.11.025 ", "10.1016/j.cell.2019.11.025"),
            ("10.1016/j.cell.2019.11.025. ", "10.1016/j.cell.2019.11.025"),
            ("DOI 10.1182/ bloodadvances.2017015073", "10.1182/bloodadvances.2017015073"),
            ("DOI:10.1167/iovs.15-18117", "10.1167/iovs.15-18117"),
            ("DOI: 10.1002/biot.201200199", "10.1002/biot.201200199"),
            ("DOI: 10.1111/j.1440-1827.1995.tb03518.x.", "10.1111/j.1440-1827.1995.tb03518.x"),
            ("https://doi.org/10.1101/2021.01.02.425073", "10.1101/2021.01.02.425073"),
            ("", ""),
        ]

        for doi, expected in test_cases:
            with self.subTest(doi=doi):
                self.assertEqual(clean_doi(doi), expected)
