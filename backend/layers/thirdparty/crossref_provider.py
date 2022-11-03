class CrossrefProviderInterface:
    def fetch_metadata(self, doi: str) -> dict:
        pass

    def fetch_preprint_published_doi(self, doi):
        pass
