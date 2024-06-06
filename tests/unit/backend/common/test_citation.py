from backend.common.citation import format_citation_crossref, format_citation_dp


def test__format_citation(self):
    message = {
        "author": [{"family": "Doe", "given": "John"}],
        "container-title": ["Test Journal"],
        "created": {"date-parts": [[2022]]},
    }
    result = format_citation_crossref(message)
    self.assertEqual(result, "Doe, John et al. (2022) Test Journal")

    message = {
        "authors": [
            {"family": "Gabitto", "given": "Mariano I."},
            {"family": "Travaglini", "given": "Kyle J."},
            {"family": "Rachleff", "given": "Victoria M."},
            {"family": "Kaplan", "given": "Eitan S."},
        ],
        "is_preprint": True,
        "journal": "bioRxiv",
        "published_at": 1683590400.0,
        "published_day": 9,
        "published_month": 5,
        "published_year": 2023,
    }
    result = format_citation_dp(message)
    self.assertEqual(result, "Gabitto, Mariano I. et al. (2023) bioRxiv")
