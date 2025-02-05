import jsonschema
import pytest

from backend.layers.common.ingestion_manifest import get_schema, get_validator, to_manifest


def test_validate_schema():
    jsonschema.Draft202012Validator.check_schema(get_schema())


@pytest.fixture(scope="module")
def ingestion_manifest_validator():
    return get_validator()


@pytest.mark.parametrize(
    "manifest",
    [
        {"anndata": "https://example.com/dataset.h5ad"},
        {"anndata": "s3://bucket/dataset.h5ad"},
        {"anndata": "https://example.com/dataset.h5ad", "atac_seq_fragments": "https://example.com/fragments.tsv.gz"},
        {"anndata": "https://example.com/dataset.h5ad", "atac_seq_fragments": "s3://bucket/fragments.tsv.gz"},
    ],
)
def test_anndata_validation_success(manifest, ingestion_manifest_validator):
    ingestion_manifest_validator.validate(manifest)


@pytest.mark.parametrize(
    "manifest",
    [{"atac_seq_fragments": "https://example.com/fragments.tsv.gz"}, {"anndata": 1234}],
)
def test_anndata_validation_failure(manifest, ingestion_manifest_validator):
    with pytest.raises(jsonschema.ValidationError):
        ingestion_manifest_validator.validate(manifest)


@pytest.mark.parametrize(
    "anndata,atac_seq_fragments",
    [
        ("https://example.com/dataset.h5ad", None),
        ("https://example.com/dataset.h5ad", "https://example.com/fragments.tsv.gz"),
    ],
)
def test_to_manifest(anndata, atac_seq_fragments, ingestion_manifest_validator):
    manifest = to_manifest(anndata, atac_seq_fragments)
    ingestion_manifest_validator.validate(manifest)
