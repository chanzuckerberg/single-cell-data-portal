from unittest import mock

from backend.layers.common.entities import CollectionVersionId


class TestCollectionMigrate:
    def test_can_open_revision_true(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        published = collections["published"][0]
        collection_version_id = CollectionVersionId()
        schema_migrate.business_logic.create_collection_version.return_value = mock.Mock(
            version_id=collection_version_id
        )
        datasets = [
            {"dataset_id": dataset.dataset_id.id, "dataset_version_id": dataset.version_id.id}
            for dataset in published.datasets
        ]
        response = schema_migrate.collection_migrate(published.collection_id.id, published.version_id.id, True)
        for i in range(len(response)):
            assert response[i]["collection_id"] == collection_version_id.id
            response[i].pop("collection_id")
            assert response[i] == datasets[i]

    def test_can_open_revision_false(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        private = collections["private"][0]
        datasets = [
            {"dataset_id": dataset.dataset_id.id, "dataset_version_id": dataset.version_id.id}
            for dataset in private.datasets
        ]
        response = schema_migrate.collection_migrate(private.collection_id.id, private.version_id.id, False)
        for i in range(len(response)):
            assert response[i]["collection_id"] == private.version_id.id
            response[i].pop("collection_id")
            assert response[i] == datasets[i]
