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
        dataset_version_id = private.datasets[0].version_id.id
        response = schema_migrate.dataset_migrate(
            private.version_id.id, private.collection_id.id, private.datasets[0].dataset_id.id, dataset_version_id
        )
        assert response["collection_version_id"] == private.version_id.id
        assert response["dataset_version_id"] == new_dataset_version_id.id
        assert dataset_version_id != new_dataset_version_id.id
        assert response["uri"] == f"s3://artifact-bucket/{dataset_version_id}/migrated.h5ad"
        assert response["sfn_name"].startswith("migrate_")
        assert new_dataset_version_id.id in response["sfn_name"]
        schema_migrate.schema_validator.migrate.assert_called_once_with(
            "previous_schema.h5ad", "migrated.h5ad", private.collection_id.id, private.datasets[0].dataset_id.id
        )

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
        dataset_version_id = published.datasets[0].version_id.id
        response = schema_migrate.dataset_migrate(
            published.version_id.id, published.collection_id.id, published.datasets[0].dataset_id.id, dataset_version_id
        )
        assert response["collection_version_id"] == published.version_id.id
        assert response["dataset_version_id"] == new_dataset_version_id.id
        assert dataset_version_id != new_dataset_version_id.id
        assert response["uri"] == f"s3://artifact-bucket/{dataset_version_id}/migrated.h5ad"
        assert response["sfn_name"].startswith("migrate_")
        assert new_dataset_version_id.id in response["sfn_name"]
        schema_migrate.schema_validator.migrate.assert_called_once_with(
            "previous_schema.h5ad", "migrated.h5ad", published.collection_id.id, published.datasets[0].dataset_id.id
        )

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
        dataset_version_id = revision.datasets[0].version_id.id
        response = schema_migrate.dataset_migrate(
            revision.version_id.id, revision.collection_id.id, revision.datasets[0].dataset_id.id, dataset_version_id
        )
        assert response["collection_version_id"] == revision.version_id.id
        assert response["dataset_version_id"] == new_dataset_version_id.id
        assert dataset_version_id != new_dataset_version_id.id
        assert response["uri"] == f"s3://artifact-bucket/{dataset_version_id}/migrated.h5ad"
        assert response["sfn_name"].startswith("migrate_")
        assert new_dataset_version_id.id in response["sfn_name"]
        schema_migrate.schema_validator.migrate.assert_called_once_with(
            "previous_schema.h5ad", "migrated.h5ad", revision.collection_id.id, revision.datasets[0].dataset_id.id
        )
