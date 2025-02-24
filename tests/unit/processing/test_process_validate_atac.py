from typing import Tuple
from unittest.mock import Mock

import pytest

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.common.entities import (
    CollectionVersion,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetId,
    DatasetVersionId,
)
from backend.layers.common.ingestion_manifest import IngestionManifest
from backend.layers.processing.process_validate_atac import ProcessValidateATAC
from tests.unit.processing.base_processing_test import BaseProcessingTest

fragment_uri_fmt = "https://datasets.cellxgene.cziscience.com/{artifact_id}.tsv.bgz"


@pytest.fixture
def setup():
    base_test = BaseProcessingTest()
    base_test.setUpClass()
    base_test.setUp()
    base_test.schema_validator.check_anndata_requires_fragment = Mock(return_value=False)
    base_test.schema_validator.validate_atac = Mock(return_value=None)
    return base_test


@pytest.fixture
def unpublished_collection(setup) -> CollectionVersion:
    return setup.generate_unpublished_collection()


@pytest.fixture
def collection_revision(setup) -> CollectionVersion:
    return setup.generate_collection_revision()


@pytest.fixture
def unpublished_dataset(unpublished_collection, setup) -> Tuple[DatasetVersionId, DatasetId]:
    new_dataset_version = setup.database_provider.create_canonical_dataset(unpublished_collection.version_id)
    setup.database_provider.add_dataset_to_collection_version_mapping(
        unpublished_collection.version_id, new_dataset_version.version_id
    )
    return new_dataset_version.version_id, new_dataset_version.dataset_id


@pytest.fixture
def dataset_with_fragment(
    unpublished_collection,
    setup,
) -> Tuple[DatasetVersionId, DatasetId]:
    new_dataset_version = setup.database_provider.create_canonical_dataset(unpublished_collection.version_id)
    setup.database_provider.add_dataset_to_collection_version_mapping(
        unpublished_collection.version_id, new_dataset_version.version_id
    )
    artifact_id = setup.database_provider.add_dataset_artifact(
        new_dataset_version.version_id, DatasetArtifactType.ATAC_FRAGMENT, "dummy"
    )
    setup.database_provider.update_dataset_artifact(artifact_id, fragment_uri_fmt.format(artifact_id=artifact_id.id))
    artifact_id = setup.database_provider.add_dataset_artifact(
        new_dataset_version.version_id, DatasetArtifactType.ATAC_INDEX, "dummy"
    )
    setup.database_provider.update_dataset_artifact(
        artifact_id, fragment_uri_fmt.format(artifact_id=artifact_id.id) + ".tbi"
    )
    return new_dataset_version.version_id, new_dataset_version.dataset_id


@pytest.fixture
def new_fragment_uri() -> str:
    return "https://www.dropbox.com/s/fake_location/test.tsv.bgz?dl=0"


@pytest.fixture
def anndata_uri() -> str:
    return "s3://fake_bucket_name/fake_key.h5ad"


def existing_fragment_uri() -> str:
    return "https://dataset_base_uri/deadbeef-36da-4643-b3d5-ee20853084ba.tsv.bgz"


@pytest.fixture
def process_validate_atac(setup):
    proc = ProcessValidateATAC(setup.business_logic, setup.uri_provider, setup.s3_provider, setup.schema_validator)
    proc.download_from_source_uri = Mock(
        side_effect=[CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME, CorporaConstants.ORIGINAL_ATAC_FRAGMENT_FILENAME]
    )
    return proc


@pytest.mark.parametrize(
    "anndata_uri",
    [
        "s3://fake_bucket_name/fake_key.h5ad",  # existing anndata
        "https://www.dropbox.com/s/fake_location/test.h5ad?dl=0",  # new anndata
    ],
)
def test_new_fragment(
    new_fragment_uri, anndata_uri, unpublished_collection, unpublished_dataset, process_validate_atac, setup
):
    """validation will succeed, status will be updated, and fragment artifacts will be uploaded"""
    dataset_version_id, _ = unpublished_dataset
    process_validate_atac.process(
        unpublished_collection.version_id,
        dataset_version_id,
        IngestionManifest(anndata=anndata_uri, atac_fragment=new_fragment_uri),
        "fake_bucket_name",
    )
    status = setup.business_logic.get_dataset_status(dataset_version_id)
    assert status.atac_status == DatasetConversionStatus.UPLOADED

    artifacts = list(setup.business_logic.get_dataset_artifacts(dataset_version_id))
    assert len(artifacts) == 2

    atac_fragment_artifact = [a for a in artifacts if a.type == DatasetArtifactType.ATAC_FRAGMENT][0]
    assert setup.s3_provider.file_exists("fake_bucket_name", atac_fragment_artifact.uri.split("/")[-1])
    assert str(atac_fragment_artifact.id) in atac_fragment_artifact.uri

    atac_frag_index_artifact = [a for a in artifacts if a.type == DatasetArtifactType.ATAC_INDEX][0]
    assert setup.s3_provider.file_exists("fake_bucket_name", atac_frag_index_artifact.uri.split("/")[-1])
    assert str(atac_fragment_artifact.id) in atac_frag_index_artifact.uri


def test_old_fragment():
    """A published fragment is used in the manifest, this will pass validation, the artifac will be copied to the new dataset version."""
    pass


def test_no_fragment_and_required():
    """A manifest is provided without a fragment, and the annotation requires one. This will fail validation."""
    dataset_version_id, _ = unpublished_dataset
    process_validate_atac.process(
        unpublished_collection.version_id,
        dataset_version_id,
        IngestionManifest(anndata=anndata_uri, atac_fragment=new_fragment_uri),
        "fake_bucket_name",
    )
    status = setup.business_logic.get_dataset_status(dataset_version_id)
    assert status.atac_status == DatasetConversionStatus.UPLOADED

    artifacts = list(setup.business_logic.get_dataset_artifacts(dataset_version_id))
    assert len(artifacts) == 2

    atac_fragment_artifact = [a for a in artifacts if a.type == DatasetArtifactType.ATAC_FRAGMENT][0]
    assert setup.s3_provider.file_exists("fake_bucket_name", atac_fragment_artifact.uri.split("/")[-1])
    assert str(atac_fragment_artifact.id) in atac_fragment_artifact.uri

    atac_frag_index_artifact = [a for a in artifacts if a.type == DatasetArtifactType.ATAC_INDEX][0]
    assert setup.s3_provider.file_exists("fake_bucket_name", atac_frag_index_artifact.uri.split("/")[-1])
    assert str(atac_fragment_artifact.id) in atac_frag_index_artifact.uri


def test_no_fragment_and_not_required():
    """A manifest is provided without a fragment, and the annotation does not require one. This will pass validation."""
    pass


def test_no_fragment_and_optional():
    """A manifest is provided without a fragment, and the annotation has an optional fragment. This will pass validation."""
    pass


def test_existing_dataset_with_fragment_removed():
    """Updating a dataset to remove the optional fragment."""
    pass


def test_existing_dataset_with_fragment_addeed():
    """Updating a dataset to add the optional fragment."""
    pass
