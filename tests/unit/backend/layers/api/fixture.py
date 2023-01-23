def generate_mock_publisher_metadata(journal_override=None):
    return {
        "authors": [{"given": "John", "family": "Doe"}, {"given": "Jane", "family": "Doe"}],
        "published_year": 2021,
        "published_month": 11,
        "published_day": 10,
        "published_at": 1636520400.0,
        "journal": journal_override or "Nature",
        "is_preprint": False,
    }
