# ruff: noqa
import os
from unittest import mock

from backend.layers.common.entities import DatasetArtifact
from tests.unit.schema_migration.pytest_fixtures import (
    private,
    published_collection,
    revision,
    schema_migrate_and_collections,
)


class TestDatasetMigrate:
    @mock.patch.dict(os.environ, {"UPLOAD_BUCKET": "upload_bucket"})
    def test_dataset_migrate(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        private = collections["private"][0]
        schema_migrate.business_logic.s3_provider.parse_s3_uri.return_value = ("fake-bucket", "object_key.h5ad")
        schema_migrate.business_logic.get_dataset_artifacts.return_value = [
            DatasetArtifact(id=None, type="raw_h5ad", uri="s3://fake-bucket/object_key.h5ad")
        ]
        with mock.patch("backend.schema_migration.migrate.cellxgene_schema"):
            dataset_version_id = private.datasets[0].version_id.id
            response = schema_migrate.dataset_migrate(
                private.collection_id.id, private.datasets[0].dataset_id.id, dataset_version_id
            )
            assert response["collection_id"] == private.collection_id.id
            assert response["dataset_version_id"] == dataset_version_id
            assert response["url"] == f"s3://upload_bucket/{dataset_version_id}/migrated.h5ad"
