from unittest import mock

from backend.layers.common.entities import CollectionVersionId


class TestCollectionMigrate:
    def test_can_publish_true(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        published = collections["published"][0]
        collection_version_id = CollectionVersionId()
        schema_migrate.business_logic.create_collection_version.return_value = mock.Mock(
            version_id=collection_version_id
        )
        datasets = [
            {"can_publish": True, "dataset_id": dataset.dataset_id.id, "dataset_version_id": dataset.version_id.id}
            for dataset in published.datasets
        ]
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, True)
        assert response["collection_version_id"] == collection_version_id.id
        assert "no_datasets" not in response
        actual_datasets = response["datasets"]
        for i in range(len(actual_datasets)):
            assert actual_datasets[i]["collection_version_id"] == collection_version_id.id
            actual_datasets[i].pop("collection_version_id")
            assert actual_datasets[i] == datasets[i]

    def test_can_publish_false(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        private = collections["private"][0]
        datasets = [
            {"can_publish": False, "dataset_id": dataset.dataset_id.id, "dataset_version_id": dataset.version_id.id}
            for dataset in private.datasets
        ]
        response = schema_migrate.collection_migrate(private.collection_id.id, private.version_id.id, False)
        assert response["collection_version_id"] == private.version_id.id
        assert "no_datasets" not in response
        actual_datasets = response["datasets"]
        for i in range(len(actual_datasets)):
            assert actual_datasets[i]["collection_version_id"] == private.version_id.id
            actual_datasets[i].pop("collection_version_id")
            assert actual_datasets[i] == datasets[i]

    def test_can_publish_false_and_no_datasets(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        published = collections["published"][0]
        published.datasets = []
        schema_migrate.business_logic.create_collection_version.return_value = mock.Mock(
            version_id=CollectionVersionId()
        )
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, False)
        assert response["collection_version_id"] == published.version_id.id
        assert not response["datasets"]
        assert response["no_datasets"] is True
