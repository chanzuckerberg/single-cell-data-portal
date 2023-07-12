from unittest import mock

from backend.layers.common.entities import DatasetArtifact, DatasetVersionId


class TestDatasetMigrate:
    def test_dataset_migrate_private(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        private = collections["private"][0]
        schema_migrate.business_logic.s3_provider.parse_s3_uri.return_value = ("fake-bucket", "object_key.h5ad")
        schema_migrate.business_logic.get_dataset_artifacts.return_value = [
            DatasetArtifact(id=None, type="raw_h5ad", uri="s3://fake-bucket/object_key.h5ad")
        ]
        new_dataset_version_id = DatasetVersionId()
        schema_migrate.business_logic.ingest_dataset.return_value = (
            new_dataset_version_id,
            private.datasets[0].dataset_id.id,
        )
        with mock.patch("backend.layers.processing.schema_migration.migrate"):
            dataset_version_id = private.datasets[0].version_id.id
            response = schema_migrate.dataset_migrate(
                private.collection_id.id, private.datasets[0].dataset_id.id, dataset_version_id
            )
            assert response["collection_version_id"] == private.collection_id.id
            assert response["dataset_version_id"] == new_dataset_version_id.id
            assert dataset_version_id != new_dataset_version_id.id
            assert response["url"] == f"s3://artifact-bucket/{dataset_version_id}/migrated.h5ad"

    def test_dataset_migrate_published(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        published = collections["published"][0]
        schema_migrate.business_logic.s3_provider.parse_s3_uri.return_value = ("fake-bucket", "object_key.h5ad")
        schema_migrate.business_logic.get_dataset_artifacts.return_value = [
            DatasetArtifact(id=None, type="raw_h5ad", uri="s3://fake-bucket/object_key.h5ad")
        ]
        new_dataset_version_id = DatasetVersionId()
        schema_migrate.business_logic.ingest_dataset.return_value = (
            new_dataset_version_id,
            published.datasets[0].dataset_id.id,
        )
        with mock.patch("backend.layers.processing.schema_migration.migrate"):
            dataset_version_id = published.datasets[0].version_id.id
            response = schema_migrate.dataset_migrate(
                published.collection_id.id, published.datasets[0].dataset_id.id, dataset_version_id
            )
            assert response["collection_version_id"] == published.collection_id.id
            assert response["dataset_version_id"] == new_dataset_version_id.id
            assert dataset_version_id != new_dataset_version_id.id
            assert response["url"] == f"s3://artifact-bucket/{dataset_version_id}/migrated.h5ad"

    def test_dataset_migrate_revision(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        revision = collections["revision"][0]
        schema_migrate.business_logic.s3_provider.parse_s3_uri.return_value = ("fake-bucket", "object_key.h5ad")
        schema_migrate.business_logic.get_dataset_artifacts.return_value = [
            DatasetArtifact(id=None, type="raw_h5ad", uri="s3://fake-bucket/object_key.h5ad")
        ]
        new_dataset_version_id = DatasetVersionId()
        schema_migrate.business_logic.ingest_dataset.return_value = (
            new_dataset_version_id,
            revision.datasets[0].dataset_id.id,
        )
        with mock.patch("backend.layers.processing.schema_migration.migrate"):
            dataset_version_id = revision.datasets[0].version_id.id
            response = schema_migrate.dataset_migrate(
                revision.collection_id.id, revision.datasets[0].dataset_id.id, dataset_version_id
            )
            assert response["collection_version_id"] == revision.collection_id.id
            assert response["dataset_version_id"] == new_dataset_version_id.id
            assert dataset_version_id != new_dataset_version_id.id
            assert response["url"] == f"s3://artifact-bucket/{dataset_version_id}/migrated.h5ad"
