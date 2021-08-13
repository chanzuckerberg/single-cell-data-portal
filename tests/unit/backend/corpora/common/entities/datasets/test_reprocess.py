from tests.unit.backend.corpora.common.entities.datasets import TestDataset


class TestReprocessDataset(TestDataset):
    def test_reprocess(self):
        def _verify(_dataset):
            previous_dataset_id = dataset.id
            previous_revision = dataset.revision
            dataset.reprocess()
            self.assertEqual(previous_dataset_id, dataset.id)
            self.assertLess(previous_revision, dataset.revision)
            # self.assertIsNone(dataset.rev)

        with self.subTest("published_dataset"):
            dataset = self.generate_dataset_with_s3_resources(
                self.session, collection_visibility="PRIVATE", published=True
            )
            s3_objects = self.get_s3_object_paths_from_dataset(dataset)
            _verify(dataset)
            for s3_objects in s3_objects:
                self.assertS3FileExists(*s3_objects)

        with self.subTest("unpublished_dataset"):
            dataset = self.generate_dataset_with_s3_resources(
                self.session, collection_visibility="PRIVATE", published=False
            )
            s3_objects = self.get_s3_object_paths_from_dataset(dataset)
            _verify(dataset)
            for s3_objects in s3_objects:
                self.assertS3FileDoesNotExist(*s3_objects)
