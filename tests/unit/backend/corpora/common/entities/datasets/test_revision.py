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
                    if key in ("id", "dataset_id"):
                        self.assertNotEqual(rev_artifacts[i][key], original_artifacts[i][key])
                    else:
                        self.assertEqual(rev_artifacts[i][key], original_artifacts[i][key])

        with self.subTest("deployment is correctly created and points to correct s3 uri "):
            rev_deployment = rev_dataset.pop("deployment_directories")[0]
            original_deployment = dataset.pop("deployment_directories")[0]
            for key in rev_deployment.keys():
                if key in ("dataset_id", "id"):
                    self.assertNotEqual(rev_deployment[key], original_deployment[key])
                else:
                    self.assertEqual(rev_deployment[key], original_deployment[key])

        with self.subTest("Test processing status copied over"):
            rev_processing_status = rev_dataset.pop("processing_status")
            original_processing_status = dataset.pop("processing_status")
            for key in rev_processing_status.keys():
                if key in ("dataset_id", "id"):
                    self.assertNotEqual(rev_processing_status[key], original_processing_status[key])
                else:
                    self.assertEqual(rev_processing_status[key], original_processing_status[key])

        with self.subTest("revision count is correctly incremented"):
            self.assertEqual(rev_dataset.pop("revision"), dataset.pop("revision") + 1)
        with self.subTest("revision points at a different collection"):
            revision_collection = rev_dataset.pop("collection")
            dataset_1_collection = dataset.pop("collection")
            self.assertNotEqual(revision_collection, dataset_1_collection)
        with self.subTest("metadata of revised matches original"):
            for key in rev_dataset.keys():
                if key in ("published", "revision", "original_id", "id", "collection_visibility"):
                    self.assertNotEqual(rev_dataset[key], dataset[key])
                else:
                    self.assertEqual(rev_dataset[key], dataset[key])
