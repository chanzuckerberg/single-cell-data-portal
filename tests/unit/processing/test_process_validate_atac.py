from typing import Tuple
from unittest.mock import Mock

import pytest

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.common.entities import (
    CollectionVersion,
    CollectionVersionWithDatasets,
    DatasetArtifact,
    DatasetArtifactId,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetId,
    DatasetStatusKey,
    DatasetValidationStatus,
    DatasetVersionId,
)
from backend.layers.common.ingestion_manifest import IngestionManifest
from backend.layers.processing.exceptions import ConversionFailed, ValidationAtacFailed
from backend.layers.processing.process_validate_atac import ProcessValidateATAC
from tests.unit.processing.base_processing_test import BaseProcessingTest

fragment_uri_fmt = "http://domain/{artifact_id}-fragment.tsv.bgz"


@pytest.fixture
def setup():
    base_test = BaseProcessingTest()
    base_test.setUpClass()
    base_test.setUp()
    base_test.schema_validator.check_anndata_requires_fragment = Mock(return_value=False)
    base_test.schema_validator.validate_atac = Mock(return_value=([], "fragment.tsv.bgz", "fragment.tsv.bgz.tbi"))
    return base_test


@pytest.fixture
def migration_set(monkeypatch):
    monkeypatch.setenv("MIGRATION", "true")


@pytest.fixture
def unpublished_collection(setup) -> CollectionVersion:
    return setup.generate_unpublished_collection()


@pytest.fixture
def unpublished_dataset(unpublished_collection, setup) -> Tuple[DatasetVersionId, DatasetId]:
    new_dataset_version = setup.database_provider.create_canonical_dataset(unpublished_collection.version_id)
    setup.database_provider.add_dataset_to_collection_version_mapping(
        unpublished_collection.version_id, new_dataset_version.version_id
    )
    return new_dataset_version.version_id, new_dataset_version.dataset_id


@pytest.fixture
def process_validate_atac(setup):
    proc = ProcessValidateATAC(setup.business_logic, setup.uri_provider, setup.s3_provider, setup.schema_validator)
    proc.download_from_source_uri = Mock(
        side_effect=[CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME, CorporaConstants.ORIGINAL_ATAC_FRAGMENT_FILENAME]
    )
    return proc


@pytest.fixture
def collection_revision_with_fragment(
    setup, unpublished_collection, process_validate_atac
) -> CollectionVersionWithDatasets:
    collection = setup.generate_published_collection(add_datasets=2)
    fragment_dataset = collection.datasets[0]
    artifact_id = process_validate_atac.create_atac_artifact(
        "anything",
        DatasetArtifactType.ATAC_FRAGMENT,
        fragment_dataset.version_id,
        "datasets",
    )
    process_validate_atac.create_atac_artifact(
        "anything",
        DatasetArtifactType.ATAC_INDEX,
        fragment_dataset.version_id,
        "datasets",
        artifact_id,
    )
    revision = setup.business_logic.create_collection_version(collection.collection_id)
    return setup.business_logic.get_collection_version(revision.version_id)


@pytest.fixture
def new_fragment_uri() -> str:
    return "https://www.dropbox.com/s/fake_location/test.tsv.bgz?dl=0"


@pytest.fixture
def anndata_uri() -> str:
    return "s3://fake_bucket_name/fake_key.h5ad"


@pytest.fixture
def manifest_with_fragment(anndata_uri, new_fragment_uri) -> IngestionManifest:
    return IngestionManifest(anndata=anndata_uri, atac_fragment=new_fragment_uri)


@pytest.fixture
def manifest_without_fragment(anndata_uri) -> IngestionManifest:
    return IngestionManifest(anndata=anndata_uri)


class TestProcessValidateAtac:
    """These tests assume that the anndata is atac, and a fragment is provided."""

    # collection revision
    ## dataset revised with optional fragment, fragment is added
    ## dataset revised without optional fragment, fragment is removed, still exists on old dataset version.
    ## dataset with required fragment and revised anndata, new anndata and same fragment
    ## dataset with revise required fragment, new fragment, old fragment still exists
    ## dataset with deleted fragment, fragment is removed

    def assert_old_fragment_replaced(self, artifacts, old_artifact_id, old_artifact_index_id, setup):
        atac_frag_index_artifact = [a for a in artifacts if a.type == DatasetArtifactType.ATAC_INDEX][0]
        assert setup.s3_provider.file_exists("datasets", atac_frag_index_artifact.uri.split("/")[-1])
        assert str(atac_frag_index_artifact.id) != str(old_artifact_index_id.id)

        atac_fragment_artifact = [a for a in artifacts if a.type == DatasetArtifactType.ATAC_FRAGMENT][0]
        assert setup.s3_provider.file_exists("datasets", atac_fragment_artifact.uri.split("/")[-1])
        assert str(atac_fragment_artifact.id) != str(old_artifact_id.id)

    def assert_new_fragment_added(self, artifacts, setup):
        atac_fragment_artifact = [a for a in artifacts if a.type == DatasetArtifactType.ATAC_FRAGMENT][0]
        assert setup.s3_provider.file_exists("datasets", atac_fragment_artifact.uri.split("/")[-1])
        assert atac_fragment_artifact.uri == f"s3://datasets/{atac_fragment_artifact.id}-fragment.tsv.bgz"

        atac_frag_index_artifact = [a for a in artifacts if a.type == DatasetArtifactType.ATAC_INDEX][0]
        assert setup.s3_provider.file_exists("datasets", atac_frag_index_artifact.uri.split("/")[-1])
        assert atac_frag_index_artifact.uri == f"s3://datasets/{atac_fragment_artifact.id}-fragment.tsv.bgz.tbi"

    def assert_artifacts_uploaded(self, setup, dataset_version_id) -> list[DatasetArtifact]:
        status = setup.business_logic.get_dataset_status(dataset_version_id)
        assert status.atac_status == DatasetConversionStatus.UPLOADED

        artifacts = setup.business_logic.get_dataset_version(dataset_version_id).artifacts
        assert len(artifacts) == 2

        return artifacts

    @pytest.mark.parametrize(
        "anndata_uri",
        [
            "s3://fake_bucket_name/fake_key.h5ad",  # existing anndata
            "https://www.dropbox.com/s/fake_location/test.h5ad?dl=0",  # new anndata
        ],
    )
    def test_new_fragment(
        self, manifest_with_fragment, unpublished_collection, unpublished_dataset, process_validate_atac, setup
    ):
        """validation will succeed, status will be updated, and fragment artifacts will be uploaded

        This covers cases where the collection is unpublished, and the datset is new.
        It also covers cases where the fragment is optional and present or required and present.
        """
        # Arrange
        dataset_version_id, _ = unpublished_dataset

        # Act
        process_validate_atac.process(
            unpublished_collection.version_id,
            dataset_version_id,
            manifest_with_fragment,
            "datasets",
        )

        # Assert
        artifacts = self.assert_artifacts_uploaded(setup, dataset_version_id)
        self.assert_new_fragment_added(artifacts, setup)

    def test_old_fragment(self, anndata_uri, collection_revision_with_fragment, process_validate_atac, setup):
        """A published fragment is used in the manifest, this will pass validation, the artifact will be copied to the
        new dataset version."""
        # Arrange
        process_validate_atac.hash_file = Mock(
            return_value="fake_hash"
        )  # mock the hash_file method to return the same value
        dataset = collection_revision_with_fragment.datasets[0]
        new_dataset_version = setup.database_provider.replace_dataset_in_collection_version(
            collection_revision_with_fragment.version_id, dataset.version_id
        )
        old_fragment_index_artifact_id = setup.get_artifact_type_from_dataset(
            dataset, DatasetArtifactType.ATAC_INDEX
        ).id
        old_fragment_artifact_id = setup.get_artifact_type_from_dataset(dataset, DatasetArtifactType.ATAC_FRAGMENT).id
        old_fragment_uri = fragment_uri_fmt.format(artifact_id=old_fragment_artifact_id)

        # Act
        process_validate_atac.process(
            collection_revision_with_fragment.version_id,
            new_dataset_version.version_id,
            IngestionManifest(anndata=anndata_uri, atac_fragment=old_fragment_uri),
            "datasets",
        )

        # Assert
        artifacts = self.assert_artifacts_uploaded(setup, new_dataset_version.version_id)

        atac_frag_index_artifact = [a for a in artifacts if a.type == DatasetArtifactType.ATAC_INDEX][0]
        assert setup.s3_provider.file_exists("datasets", atac_frag_index_artifact.uri.split("/")[-1])
        assert str(atac_frag_index_artifact.id) == str(old_fragment_index_artifact_id.id)

        atac_fragment_artifact = [a for a in artifacts if a.type == DatasetArtifactType.ATAC_FRAGMENT][0]
        assert setup.s3_provider.file_exists("datasets", atac_fragment_artifact.uri.split("/")[-1])
        assert str(atac_fragment_artifact.id) == str(old_fragment_artifact_id)

    def test_old_fragment_replaced_because_hash_difference(
        self, anndata_uri, collection_revision_with_fragment, process_validate_atac, setup, migration_set
    ):
        """A published fragment is used in the manifest. This will pass validation, but the hash of the new file is
        different, so a new artifact will be added to the dataset version."""
        # Arrange
        process_validate_atac.hash_file = Mock(
            side_effect=["abcd", "efgh"]
        )  # mock the hash_file method to return different values
        dataset = collection_revision_with_fragment.datasets[0]
        new_dataset_version = setup.database_provider.replace_dataset_in_collection_version(
            collection_revision_with_fragment.version_id, dataset.version_id
        )
        old_fragment_index_artifact_id = setup.get_artifact_type_from_dataset(
            dataset, DatasetArtifactType.ATAC_INDEX
        ).id
        old_fragment_artifact_id = setup.get_artifact_type_from_dataset(dataset, DatasetArtifactType.ATAC_FRAGMENT).id
        old_fragment_uri = fragment_uri_fmt.format(artifact_id=old_fragment_artifact_id)

        # Act
        process_validate_atac.process(
            collection_revision_with_fragment.version_id,
            new_dataset_version.version_id,
            IngestionManifest(anndata=anndata_uri, atac_fragment=old_fragment_uri),
            "datasets",
        )

        # Assert
        artifacts = self.assert_artifacts_uploaded(setup, new_dataset_version.version_id)
        self.assert_old_fragment_replaced(artifacts, old_fragment_artifact_id, old_fragment_index_artifact_id, setup)

    def test_replace_existing_fragment(
        self, collection_revision_with_fragment, process_validate_atac, setup, manifest_with_fragment
    ):
        # Arrange
        dataset = collection_revision_with_fragment.datasets[0]
        new_dataset_version = setup.database_provider.replace_dataset_in_collection_version(
            collection_revision_with_fragment.version_id, dataset.version_id
        )
        old_fragment_index_artifact_id = setup.get_artifact_type_from_dataset(
            dataset, DatasetArtifactType.ATAC_INDEX
        ).id
        old_fragment_artifact_id = setup.get_artifact_type_from_dataset(dataset, DatasetArtifactType.ATAC_FRAGMENT).id

        # Act
        process_validate_atac.process(
            collection_revision_with_fragment.version_id,
            new_dataset_version.version_id,
            manifest_with_fragment,
            "datasets",
        )

        # Assert
        artifacts = self.assert_artifacts_uploaded(setup, new_dataset_version.version_id)
        self.assert_old_fragment_replaced(artifacts, old_fragment_artifact_id, old_fragment_index_artifact_id, setup)

    def test_existing_dataset_with_fragment_removed(
        self, collection_revision_with_fragment, process_validate_atac, setup, manifest_without_fragment
    ):
        """Updating a dataset to remove the optional fragment."""
        # Arrange
        process_validate_atac.schema_validator.check_anndata_requires_fragment = Mock(return_value=False)
        dataset_version_id = collection_revision_with_fragment.datasets[0].version_id
        new_dataset_version = setup.database_provider.replace_dataset_in_collection_version(
            collection_revision_with_fragment.version_id, dataset_version_id
        )

        # Act
        process_validate_atac.process(
            collection_revision_with_fragment.version_id,
            new_dataset_version.version_id,
            manifest_without_fragment,
            "datasets",
        )

        # Assert
        status = setup.business_logic.get_dataset_status(new_dataset_version.version_id)
        assert status.atac_status == DatasetConversionStatus.SKIPPED

        artifacts = setup.business_logic.get_dataset_version(new_dataset_version.version_id).artifacts
        assert len(artifacts) == 0

    def test_existing_dataset_with_fragment_added(
        self, collection_revision_with_fragment, process_validate_atac, setup, manifest_with_fragment
    ):
        """Updating a dataset to add the optional fragment."""
        # Arrange
        process_validate_atac.schema_validator.check_anndata_requires_fragment = Mock(return_value=False)
        dataset_version_id = collection_revision_with_fragment.datasets[0].version_id
        new_dataset_version = setup.database_provider.replace_dataset_in_collection_version(
            collection_revision_with_fragment.version_id, dataset_version_id
        )

        # Act
        process_validate_atac.process(
            collection_revision_with_fragment.version_id,
            new_dataset_version.version_id,
            manifest_with_fragment,
            "datasets",
        )

        # Assert
        artifacts = self.assert_artifacts_uploaded(setup, new_dataset_version.version_id)
        self.assert_new_fragment_added(artifacts, setup)


class TestSkipATACValidation:
    def test_not_atac_and_no_fragment(
        self, process_validate_atac, unpublished_dataset, setup, manifest_without_fragment
    ):
        """The anndata file is not ATAC, and no fragment file in the manifest, so validation should be skipped."""
        # Arrange
        dataset_version_id, _ = unpublished_dataset
        process_validate_atac.schema_validator.check_anndata_requires_fragment = Mock(side_effect=ValueError("test"))

        # Act
        assert process_validate_atac.skip_atac_validation("fake_path", manifest_without_fragment, dataset_version_id)

        # Assert
        dataset_status = setup.business_logic.get_dataset_status(dataset_version_id)
        assert setup.business_logic.get_dataset_status(dataset_version_id).atac_status == DatasetConversionStatus.NA
        assert dataset_status.validation_message == "test"

    def test_not_atac_and_fragment(self, process_validate_atac, unpublished_dataset, setup, manifest_with_fragment):
        """A manifest is provided with a fragment, and the anndata does not require one. This will fail
        validation."""
        # Arrange
        dataset_version_id, _ = unpublished_dataset
        process_validate_atac.schema_validator.check_anndata_requires_fragment = Mock(side_effect=ValueError("test"))

        # Act
        with pytest.raises(ValidationAtacFailed) as e:
            process_validate_atac.skip_atac_validation("fake_path", manifest_with_fragment, dataset_version_id)
        # Assert
        assert e.value.errors == ["test", "Fragment file not allowed for non atac anndata."]
        dataset_status = setup.business_logic.get_dataset_status(dataset_version_id)
        assert dataset_status.validation_status == DatasetValidationStatus.INVALID

    def test_optional_and_no_fragment(
        self, process_validate_atac, unpublished_dataset, setup, manifest_without_fragment
    ):
        """A manifest is provided without a fragment, and the anndata does not require one. This will pass
        validation."""
        # Arrange
        dataset_version_id, _ = unpublished_dataset

        # Act
        assert process_validate_atac.skip_atac_validation("fake_path", manifest_without_fragment, dataset_version_id)
        # Assert
        dataset_status = setup.business_logic.get_dataset_status(dataset_version_id)
        assert (
            setup.business_logic.get_dataset_status(dataset_version_id).atac_status == DatasetConversionStatus.SKIPPED
        )
        assert dataset_status.validation_message == "Fragment is optional and not present."

    def test_optional_and_fragment(self, process_validate_atac, unpublished_dataset, setup, manifest_with_fragment):
        """A manifest is provided without a fragment, and the anndata has an optional fragment. This will pass
        validation."""
        # Arrange
        dataset_version_id, _ = unpublished_dataset

        # Act
        assert not process_validate_atac.skip_atac_validation("fake_path", manifest_with_fragment, dataset_version_id)

    def test_required_and_missing_fragment(
        self, process_validate_atac, unpublished_dataset, setup, manifest_without_fragment
    ):
        """A manifest is provided without a fragment, and the anndata requires one. This will fail validation."""
        # Arrange
        dataset_version_id, _ = unpublished_dataset
        process_validate_atac.schema_validator.check_anndata_requires_fragment = Mock(return_value=True)

        # Act
        with pytest.raises(ValidationAtacFailed) as e:
            process_validate_atac.skip_atac_validation("fake_path", manifest_without_fragment, dataset_version_id)

        # Assert
        assert e.value.errors == ["Anndata requires fragment file"]
        dataset_status = setup.business_logic.get_dataset_status(dataset_version_id)
        assert dataset_status.validation_status == DatasetValidationStatus.INVALID


class TestHashFile:
    def test_hash_file(self, process_validate_atac, tmpdir):
        """Test that the hash_file method returns the correct hash."""

        # Arrange
        file_path = tmpdir.join("test.txt")
        with open(tmpdir.join("test.txt"), "w") as f:
            f.write("test")

        # Act
        assert isinstance(process_validate_atac.hash_file(file_path), str)


class TestCreateAtacArtifact:
    def test_fragment(self, process_validate_atac, unpublished_dataset, setup):
        """Test that the create_atac_artifact method creates an artifact for the fragment."""
        # Arrange
        dataset_version_id, _ = unpublished_dataset
        artifact_id = process_validate_atac.create_atac_artifact(
            "anything",
            DatasetArtifactType.ATAC_FRAGMENT,
            dataset_version_id,
            "datasets",
        )

        # Assert
        dataset = setup.business_logic.get_dataset_version(dataset_version_id)
        assert dataset.status.atac_status == DatasetConversionStatus.UPLOADED

        artifacts = dataset.artifacts
        assert len(artifacts) == 1
        assert str(artifact_id.id) == str(artifacts[0].id)
        assert artifacts[0].type == DatasetArtifactType.ATAC_FRAGMENT
        assert artifacts[0].uri == f"s3://datasets/{artifacts[0].id}-fragment.tsv.bgz"

    def test_fragment_index(self, process_validate_atac, unpublished_dataset, setup):
        """Test that the create_atac_artifact method creates an artifact for the fragment index."""
        # Arrange
        dataset_version_id, _ = unpublished_dataset
        fragment_artifact_id = DatasetArtifactId("deadbeef-36da-4643-b3d5-ee20853084ba")
        artifact_id = process_validate_atac.create_atac_artifact(
            "anything",
            DatasetArtifactType.ATAC_INDEX,
            dataset_version_id,
            "datasets",
            fragment_artifact_id=fragment_artifact_id,
        )

        # Assert
        dataset = setup.business_logic.get_dataset_version(dataset_version_id)
        assert dataset.status.atac_status == DatasetConversionStatus.UPLOADED

        artifacts = dataset.artifacts
        assert len(artifacts) == 1
        assert str(artifact_id.id) == str(artifacts[0].id)
        assert artifacts[0].type == DatasetArtifactType.ATAC_INDEX
        assert artifacts[0].uri == f"s3://datasets/{fragment_artifact_id.id}-fragment.tsv.bgz.tbi"

    def test_exception(self, process_validate_atac, unpublished_dataset, setup):
        """Test that the create_atac_artifact method raises an exception when the artifact cannot be created."""
        # Arrange
        dataset_version_id, _ = unpublished_dataset
        process_validate_atac.business_logic.add_dataset_artifact = Mock(side_effect=ValueError("test"))

        # Act
        with pytest.raises(ConversionFailed) as e:
            process_validate_atac.create_atac_artifact(
                "anything",
                DatasetArtifactType.ATAC_FRAGMENT,
                dataset_version_id,
                "datasets",
            )

        assert e.value.failed_status == DatasetStatusKey.ATAC
