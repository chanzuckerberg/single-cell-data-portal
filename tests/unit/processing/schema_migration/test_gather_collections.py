# ruff: noqa


class TestGatherCollections:
    def test_with_revision(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        published, revision = collections["revision"]
        # get_collections is called twice
        schema_migrate.business_logic.get_collections.side_effect = [[published], [revision]]
        response = schema_migrate.gather_collections()
        assert published.collection_id.id not in response["published"]
        assert revision.version_id.id in response["revision"]
        assert [
            {"dataset_id": ds.dataset_id.id, "dataset_version_id": ds.version_id.id} for ds in revision.datasets
        ] == response["revision"][revision.version_id.id]

    def test_with_published(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        published = collections["published"]
        # get_collections is called twice
        schema_migrate.business_logic.get_collections.side_effect = [published, []]
        response = schema_migrate.gather_collections()
        assert published[0].collection_id.id in response["published"]
        assert [
            {"dataset_id": ds.dataset_id.id, "dataset_version_id": ds.version_id.id} for ds in published[0].datasets
        ] == response["published"][published[0].collection_id.id]

    def test_with_private(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        private = collections["private"]
        # get_collections is called twice
        schema_migrate.business_logic.get_collections.side_effect = [[], private]
        response = schema_migrate.gather_collections()
        assert private[0].collection_id.id in response["private"]
        assert [
            {"dataset_id": ds.dataset_id.id, "dataset_version_id": ds.version_id.id} for ds in private[0].datasets
        ] == response["private"][private[0].collection_id.id]
