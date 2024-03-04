from backend.layers.processing.utils import rds_citation_from_h5ad


class TestUtils:
    def test_rds_citation_from_h5ad(self):
        # No-op
        citation_without_h5ad = "Publication: https://doi.org/10.1126/science.abl4896 Dataset Version: https://datasets.cellxgene.cziscience.com/dbd8b789-3efa-4a63-9243-90cff64f2045.rds curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5"
        assert rds_citation_from_h5ad(citation_without_h5ad) == citation_without_h5ad

        # Simple case: replace <dataset_id>.h5ad with <dataset_id>.rds
        citation_with_h5ad = "Publication: https://doi.org/10.1126/science.abl4896 Dataset Version: https://datasets.cellxgene.cziscience.com/dbd8b789-3efa-4a63-9243-90cff64f2045.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5"
        assert (
            rds_citation_from_h5ad(citation_with_h5ad)
            == "Publication: https://doi.org/10.1126/science.abl4896 Dataset Version: https://datasets.cellxgene.cziscience.com/dbd8b789-3efa-4a63-9243-90cff64f2045.rds curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5"
        )

        # Ensure we're also not replacing the .h5ad in the DOI, since it does not precede a dataset ID
        citation_with_h5ad_in_doi = "Publication: https://doi.org/10.1126/science.h5ad Dataset Version: https://datasets.cellxgene.cziscience.com/dbd8b789-3efa-4a63-9243-90cff64f2045.h5ad curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5"
        assert (
            rds_citation_from_h5ad(citation_with_h5ad_in_doi)
            == "Publication: https://doi.org/10.1126/science.h5ad Dataset Version: https://datasets.cellxgene.cziscience.com/dbd8b789-3efa-4a63-9243-90cff64f2045.rds curated and distributed by CZ CELLxGENE Discover in Collection: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5"
        )
