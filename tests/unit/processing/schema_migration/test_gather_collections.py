class TestGatherCollections:
    def test_with_revision(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        published, revision = collections["revision"]
        # get_collections is called twice
        schema_migrate.business_logic.get_collections.side_effect = [[published], [revision]]
        response = schema_migrate.gather_collections()
        assert {"can_publish": "True", "collection_id": published.collection_id.id} not in response
        assert {
            "can_publish": "False",
            "collection_id": revision.collection_id.id,
            "collection_version_id": revision.version_id.id,
        } in response
        assert len(response) == 1

    def test_with_published(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        published = collections["published"]
        # get_collections is called twice
        schema_migrate.business_logic.get_collections.side_effect = [published, []]
        response = schema_migrate.gather_collections()
        assert {
            "can_publish": "True",
            "collection_id": published[0].collection_id.id,
            "collection_version_id": published[0].version_id.id,
        } in response
        assert len(response) == 1

    def test_with_private(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        private = collections["private"]
        # get_collections is called twice
        schema_migrate.business_logic.get_collections.side_effect = [[], private]
        response = schema_migrate.gather_collections()
        assert {
            "can_publish": "False",
            "collection_id": private[0].collection_id.id,
            "collection_version_id": private[0].version_id.id,
        } in response
        assert len(response) == 1
