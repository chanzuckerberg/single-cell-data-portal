import unittest

from backend.corpora.common.entities.tiledb_data import TileDBData


class TestTileDbData(unittest.TestCase):
    def setUp(self):
        self.metadata = {
            "name": "test_collection",
            "description": "This is a test collection.",
            "owner": "Rohan Agarwal",
            "contact_name": "Emanuele Bezzi",
            "contact_email": "ebezzi@chanzuckerberg.com",
            "curator_name": "Kuni Katsuya",
            "links": [
                {"link_url": "https://doi.org/fake-article", "link_name": "Fake", "link_type": "FAKE TYPE"},
                {"link_url": "https://doi.org/test-article", "link_name": "Test", "link_type": "Test TYPE"},
            ],
            "publisher_metadata": {
                "authors": [
                    {
                        "family": "Agarwal",
                        "given": "Rohan",
                        "name": ""
                    },
                ],
                "is_preprint": True,
                "journal": "Nature",
                "published_day": 10,
                "published_month": 10,
                "published_year": 2022
            }
        }
        self.dataset_metadata = {
            "x_approximate_distribution": "COUNT",
            "x_normalization": "string",
            "assay": [
                {
                    "label": "string",
                    "ontology_term_id": "string"
                }
            ],
            "cell_count": 12,
            "cell_type": [
                {
                    "label": "string",
                    "ontology_term_id": "string"
                }
            ],
            "development_stage": [
                {
                    "label": "string",
                    "ontology_term_id": "string"
                }
            ],
            "disease": [
                {
                    "label": "string",
                    "ontology_term_id": "string"
                }
            ],
            "ethnicity": [
                {
                    "label": "string",
                    "ontology_term_id": "string"
                }
            ],
            "dataset_assets": [
                {
                    "dataset_id": "string",
                    "filename": "string",
                    "filetype": "H5AD",
                    "id": "string",
                    "s3_uri": "string"
                }
            ],
            "is_primary_data": "PRIMARY",
            "name": "string",
            "organism": [
                {
                    "label": "string",
                    "ontology_term_id": "string"
                }
            ],
            "sex": [
                {
                    "label": "string",
                    "ontology_term_id": "string"
                }
            ],
            "tissue": [
                {
                    "label": "string",
                    "ontology_term_id": "string"
                }
            ],
            "explorer_url": "https://cellxgene.com/e/lolol.cxg",
            "processing_status": {
                "processing_status": "PENDING"
            }
        }
        self.location = "tests/unit/backend/fixtures/test_tiledb/metadata"
        TileDBData.init_db(self.location)
        self.db = TileDBData(self.location)

    def tearDown(self):
        TileDBData.destroy_db(self.location)

    def create_collection(self):
        id = self.db.create_collection({
            "name": self.metadata['name'],
            "description": self.metadata['description'],
            "owner": self.metadata['owner'],
            "contact_name": self.metadata['contact_name'],
            "contact_email": self.metadata['contact_email'],
            "curator_name": self.metadata['curator_name'],
            "links": self.metadata['links'],
            "publisher_metadata": self.metadata['publisher_metadata']
        })
        return id

    def test_create_collection(self):
        id = self.create_collection()
        coll = self.db.get_collection(id)
        self.assertEqual(coll['name'], self.metadata['name'])
        self.assertDictEqual(coll['publisher_metadata'], self.metadata['publisher_metadata'])

    def test_edit_collection(self):
        # Create collection
        id = self.create_collection()
        # Edit collection name
        self.db.edit_collection(id, "name", "edited_name")
        coll = self.db.get_collection(id)
        self.assertEqual(coll['name'], "edited_name")
        # Edit publisher metadata
        new_data = {
                "authors": [
                    {
                        "family": "Agarwal",
                        "given": "Rohan",
                        "name": ""
                    },
                ],
                "is_preprint": True,
                "journal": "Nature",
                "published_day": 10,
                "published_month": 10,
                "published_year": 2022
            }
        self.db.edit_collection(id, "publisher_metadata", new_data)
        coll = self.db.get_collection(id)
        self.assertDictEqual(coll['publisher_metadata'], new_data)

    def test_add_dataset_to_collection(self):
        # Create collection
        id = self.create_collection()
        # Add dataset to collection
        dataset_id = self.db.add_dataset(id, self.dataset_metadata)
        # Get dataset/collection
        datasets = self.db.get_datasets(id)
        dataset_ids = self.db.get_attribute(id, "datasets")
        self.assertEqual(datasets[0]['cell_count'], 12)
        self.assertEqual(dataset_ids[0], dataset_id)

    def test_delete_dataset_from_collection(self):
        # Create collection
        id = self.create_collection()
        # Add dataset to collection
        dataset_id_1 = self.db.add_dataset(id, self.dataset_metadata)
        dataset_id_2 = self.db.add_dataset(id, self.dataset_metadata)
        datasets = self.db.get_datasets(id)
        self.assertEqual(len(datasets), 2)
        # Get datasets and delete one
        self.db.delete_dataset(id, dataset_id_1)
        # Get dataset/collection
        datasets = self.db.get_datasets(id)
        self.assertEqual(len(datasets), 1)
        self.assertEqual(datasets[0]['id'][0], dataset_id_2)

    def test_collection_history(self):
        # create collection
        id = self.create_collection()
        # edit collection
        self.db.edit_collection(id, "name", "edited_name")
        # read previous version
        coll = self.db.read_collection_history(id, 1)  # 1 version back
        self.assertEqual(coll['name'], "test_collection")
        # revert to a previous version
        self.db.revert_collection_history(id, 1)
        coll = self.db.get_collection(id)
        self.assertEqual(coll['name'], "test_collection")

    def test_dataset_history(self):
        # create collection
        id = self.create_collection()
        # add dataset
        did = self.db.add_dataset(id, self.dataset_metadata)
        # edit dataset
        self.db.edit_dataset(did, "cell_count", 13)
        # read previous version
        dataset = self.db.read_dataset_history(did, 1)  # 1 version back
        self.assertEqual(dataset['cell_count'], 12)
        # revert to previous version
        self.db.revert_dataset_history(did, 1)
        dataset = self.db.get_dataset(did)
        self.assertEqual(dataset['cell_count'], 12)

    def test_create_revision_on_collection(self):
        # Create a collection
        id = self.create_collection()
        # Create the revision
        rev_id = self.db.create_revision(id)
        self.assertNotEqual(rev_id, id)
        # Get the revision
        rev = self.db.get_collection(rev_id)
        self.assertEqual(rev['name'], "test_collection")
        self.assertEqual(rev['visibility'], "PRIVATE")

    def test_edit_revision(self):
        # Create a collection
        id = self.create_collection()
        # Create a revision
        rev_id = self.db.create_revision(id)
        # Edit the revision
        self.db.edit_collection(rev_id, "name", "edited_name")
        # Get the revision
        rev = self.db.get_collection(rev_id)
        self.assertEqual(rev['name'], "edited_name")

    def test_publish_revision(self):
        # Create a collection
        id = self.create_collection()
        # Create a revision
        rev_id = self.db.create_revision(id)
        self.db.edit_collection(rev_id, "name", "edited_name")
        # Publish the revision
        self.db.publish_revision(rev_id)
        # Get the new collection
        coll = self.db.get_collection(id)
        self.assertEqual(coll['name'], "edited_name")
        self.assertEqual(coll['visibility'], "PUBLIC")

    def test_delete_revision(self):
        # Create a collection
        id = self.create_collection()
        # Create a revision
        rev_id = self.db.create_revision(id)
        # Delete the revision
        self.db.delete_collection(rev_id)
        # Try getting the revision
        status = self.db.get_attribute(rev_id, "visibility")
        self.assertEqual(status, "DELETED")

    def test_get_published_collections(self):
        _ = self.create_collection()
        id2 = self.create_collection()
        self.db.publish_collection(id2)
        colls = self.db.get_published_collections()
        self.assertEqual(colls[0]['id'], id2)

    def test_get_published_datasets(self):
        id1 = self.create_collection()
        _ = self.db.add_dataset(id1, self.dataset_metadata)
        id2 = self.create_collection()
        d2 = self.db.add_dataset(id2, self.dataset_metadata)
        self.db.publish_collection(id2)
        datasets = self.db.get_published_datasets()
        self.assertEqual(datasets[0]['id'], d2)

    def test_check_collection_access(self):
        id = self.create_collection()
        self.db.edit_collection(id, "visibility", "PRIVATE")
        self.assertTrue(self.db.check_collection_access(id, "Rohan Agarwal"))
        self.assertFalse(self.db.check_collection_access(id, "Emanuele Bezzi"))
        self.db.edit_collection(id, "visibility", "PUBLIC")
        self.assertTrue(self.db.check_collection_access(id, "Emanuele Bezzi"))


if __name__ == "__main__":
    unittest.main()
