import pydantic
import pytest

from backend.layers.common.ingestion_manifest import IngestionManifest


@pytest.mark.parametrize(
    "manifest",
    [
        """{"anndata": "https://example.com/dataset.h5ad"}""",
        """{"anndata": "s3://bucket/dataset.h5ad"}""",
        """{"anndata": "https://example.com/dataset.h5ad", "atac_seq_fragments": "https://example.com/fragments.tsv.gz"}""",
        """{"anndata": "https://example.com/dataset.h5ad", "atac_seq_fragments": "s3://bucket/fragments.tsv.gz"}""",
    ],
)
def test_anndata_validation_success(manifest):
    IngestionManifest.model_validate_json(manifest)


@pytest.mark.parametrize(
    "manifest",
    ["""{"atac_seq_fragments": "https://example.com/fragments.tsv.gz"}""", """{"anndata": 1234}"""],
)
def test_anndata_validation_failure(manifest):
    with pytest.raises(pydantic.ValidationError):
        IngestionManifest.model_validate_json(manifest)
