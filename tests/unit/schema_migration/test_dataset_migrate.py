import os
from unittest import mock

from backend.layers.common.entities import DatasetArtifact, DatasetVersionId


class TestDatasetMigrate:
    @mock.patch.dict(os.environ, {"UPLOAD_BUCKET": "upload_bucket"})
    def test_dataset_migrate(self, schema_migrate_and_collections):
        schema_migrate, collections = schema_migrate_and_collections
        private = collections["private"][0]
        schema_migrate.business_logic.s3_provider.parse_s3_uri.return_value = ("fake-bucket", "object_key.h5ad")
        schema_migrate.business_logic.get_dataset_artifacts.return_value = [
            DatasetArtifact(id=None, type="raw_h5ad", uri="s3://fake-bucket/object_key.h5ad")
        ]
        dataset_version_id = DatasetVersionId().id
        schema_migrate.business_logic.ingest_dataset.return_value = (
            dataset_version_id,
            private.datasets[0].dataset_id.id,
        )
        with mock.patch("backend.schema_migration.migrate.cellxgene_schema"):
            response = schema_migrate.dataset_migrate(
                private.collection_id.id, private.datasets[0].dataset_id.id, private.datasets[0].version_id.id
            )
            assert response["new_dataset_version_id"] == dataset_version_id
