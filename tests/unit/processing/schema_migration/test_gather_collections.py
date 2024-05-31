class TestGatherCollections:
    def test_with_revision(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections

        published, revision = collections["revision"]
        # get_collections is called twice
        schema_migrate.business_logic.get_collections.side_effect = [[published], [revision]]
        _, response = schema_migrate.gather_collections(auto_publish=True)
        assert {
            "can_publish": "True",
            "collection_id": published.collection_id.id,
            "collection_version_id": published.collection_id.id,
            "execution_id": "test-execution-arn",
        } not in response
        assert {
            "can_publish": "False",
            "collection_id": revision.collection_id.id,
            "collection_version_id": revision.version_id.id,
            "execution_id": "test-execution-arn",
        } in response
        assert len(response) == 1

    def test_with_published(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        published = collections["published"]
        # get_collections is called twice
        schema_migrate.business_logic.get_collections.side_effect = [published, []]
        _, response = schema_migrate.gather_collections(auto_publish=True)
        assert {
            "can_publish": "True",
            "collection_id": published[0].collection_id.id,
            "collection_version_id": published[0].version_id.id,
            "execution_id": "test-execution-arn",
        } in response
        assert len(response) == 1

    def test_with_private(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        private = collections["private"]
        # get_collections is called twice
        schema_migrate.business_logic.get_collections.side_effect = [[], private]
        _, response = schema_migrate.gather_collections(auto_publish=True)
        assert {
            "can_publish": "False",
            "collection_id": private[0].collection_id.id,
            "collection_version_id": private[0].version_id.id,
            "execution_id": "test-execution-arn",
        } in response
        assert len(response) == 1

    def test_with_auto_publish_false(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        published_no_revision = collections["published"][0]
        private = collections["private"][0]
        published_with_revision, revision = collections["revision"]
        # get_collections is called twice
        schema_migrate.business_logic.get_collections.side_effect = [
            [published_no_revision, published_with_revision],
            [private, revision],
        ]
        _, response = schema_migrate.gather_collections(auto_publish=False)
        for obj in [
            {
                "can_publish": "False",
                "collection_id": published_no_revision.collection_id.id,
                "collection_version_id": published_no_revision.version_id.id,
                "execution_id": "test-execution-arn",
            },
            {
                "can_publish": "False",
                "collection_id": private.collection_id.id,
                "collection_version_id": private.version_id.id,
                "execution_id": "test-execution-arn",
            },
            {
                "can_publish": "False",
                "collection_id": revision.collection_id.id,
                "collection_version_id": revision.version_id.id,
                "execution_id": "test-execution-arn",
            },
        ]:
            assert obj in response
