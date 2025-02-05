import jsonschema
import pytest

from backend.layers.common.ingest_manifest import schema, to_manifest, validator


def test_validate_schema():
    jsonschema.Draft7Validator.check_schema(schema)


@pytest.mark.parametrize(
    "manifest",
    [
        {"anndata": "https://example.com/dataset.h5ad"},
        {"anndata": "s3://bucket/dataset.h5ad"},
        {"anndata": "https://example.com/dataset.h5ad", "atac_seq_fragments": "https://example.com/fragments.tsv.gz"},
        {"anndata": "https://example.com/dataset.h5ad", "atac_seq_fragments": "s3://bucket/fragments.tsv.gz"},
    ],
)
def test_anndata_validation_success(manifest):
    validator.validate(manifest)


@pytest.mark.parametrize(
    "manifest",
    [
        {"atac_seq_fragments": "https://example.com/fragments.tsv.gz"},
        {"anndata": 1234},
    ],
)
def test_anndata_validation_failure(manifest):
    with pytest.raises(jsonschema.ValidationError):
        validator.validate(manifest)


@pytest.mark.parametrize(
    "anndata,atac_seq_fragments",
    [
        ("https://example.com/dataset.h5ad", None),
        ("https://example.com/dataset.h5ad", "https://example.com/fragments.tsv.gz"),
    ],
)
def test_to_manifest(anndata, atac_seq_fragments):
    manifest = to_manifest(anndata, atac_seq_fragments)
    validator.validate(manifest)
