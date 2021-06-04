from tests.unit.backend.corpora.common.entities.datasets import TestDataset


class TestDatasetRevision(TestDataset):
    def test__create_dataset_revision(self):
        dataset = self.generate_dataset_with_s3_resources(self.session, published=True)
        rev_dataset = dataset.create_revision().to_dict()
        dataset = dataset.to_dict()

        with self.subTest("artifacts are correctly created and point to correct s3 uri"):
            rev_artifacts = rev_dataset.pop("artifacts")
            original_artifacts = dataset.pop("artifacts")
            for i in range(0, len(rev_artifacts)):
                for key in rev_artifacts[i].keys():
                    self.compare_original_and_revision(
                        original_artifacts[i], rev_artifacts[i], key, ("dataset_id", "id")
                    )

        with self.subTest("deployment is correctly created and points to correct s3 uri "):
            rev_deployment = rev_dataset.pop("explorer_url")
            original_deployment = dataset.pop("explorer_url")
            self.assertEqual(rev_deployment, original_deployment)

        with self.subTest("Test processing status copied over"):
            rev_processing_status = rev_dataset.pop("processing_status")
            original_processing_status = dataset.pop("processing_status")
            for key in rev_processing_status.keys():
                self.compare_original_and_revision(
                    original_processing_status, rev_processing_status, key, ("dataset_id", "id")
                )

        with self.subTest("revision points at a different collection"):
            revision_collection = rev_dataset.pop("collection")
            dataset_1_collection = dataset.pop("collection")
            self.assertNotEqual(revision_collection, dataset_1_collection)
        with self.subTest("metadata of revised matches original"):
            for key in rev_dataset.keys():
                self.compare_original_and_revision(
                    dataset, rev_dataset, key, ("original_id", "id", "collection_visibility")
                )

    def compare_original_and_revision(self, original, revision, key, unique_fields):
        if key in unique_fields:
            self.assertNotEqual(original[key], revision[key])
        else:
            self.assertEqual(original[key], revision[key])
