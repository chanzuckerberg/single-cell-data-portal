import unittest

from backend.corpora.common.entities.tiledb_data import TileDBData

class TestSession(unittest.TestCase):
    def setUp(self):
        self.metadata = {
            "name": "test_collection",
            "description": "This is a test collection.",
            "owner": "Rohan Agarwal",
            "contact_name": "Emanuele Bezzi",
            "contact_email": "ebezzi@chanzuckerberg.com",
            "links": ["https://nature.com/fake-article", "https://science.com/test-paper"],
        } # TODO: fill in rest of required metadata
        self.location = "/Users/ragarwal/code/single-cell-data-portal/tests/unit/backend/fixtures/test_tiledb/metadata"
        TileDBData.init_db(self.location)
        self.db = TileDBData(self.location)

    def tearDown(self):
        TileDBData.destroy_db(self.location)

    def createCollection(self):
        id = self.db.create_collection(
            self.metadata['name'], 
            self.metadata['description'], 
            self.metadata['owner'],
            self.metadata['contact_name'], 
            self.metadata['contact_email'], 
            self.metadata['links'],
            datasets=[]
        )
        return id

    def test_create_collection(self):
        id = self.createCollection()
        coll = self.db.get_collection(id)
        self.assertEqual(coll['name'][0], self.metadata['name']) # TODO: figure out why attributes are stored as arrays

    def test_edit_collection(self):
        # Create collection 
        id = self.createCollection()
        # Edit collection name
        self.db.edit_collection(id, "name", "edited_name")
        # check
        coll = self.db.get_collection(id)
        self.assertEqual(coll['name'][0], "edited_name")
        # TODO: test editing other fields

    def test_add_dataset_to_collection(self):
        # Create collection 
        id = self.createCollection()
        # Add dataset to collection
        dataset_name = "test_dataset"
        dataset_id = self.db.add_dataset(id, dataset_name)
        # Get dataset/collection
        datasets = self.db.get_datasets(id)
        dataset_ids = self.db.get_attribute(id, "datasets")
        self.assertEqual(datasets[0]['name'], dataset_name)
        self.assertEqual(dataset_ids[0], dataset_id)
        # TODO: artifacts stuff

    def test_delete_dataset_from_collection(self):
        # Create collection 
        id = self.createCollection()
        # Add dataset to collection
        dataset_name_1 = "test_1"
        dataset_name_2 = "test_2"
        self.db.add_dataset(id, dataset_name_1)
        self.db.add_dataset(id, dataset_name_2)
        # Get datasets and delete one
        datasets = self.db.get_datasets(id)
        to_del = datasets[0]['uuid'][0].decode("utf-8") 
        to_keep = datasets[1]['uuid'][0]
        self.db.delete_dataset(id, to_del)
        # Get dataset/collection
        datasets = self.db.get_datasets(id)
        self.assertEqual(datasets[0]['uuid'][0], to_keep)

    def test_collection_history(self):
        # create collection
        id = self.createCollection()
        # edit collection
        self.db.edit_collection(id, "name", "edited_name")
        # read previous version
        coll = self.db.read_collection_history(id, 1) # 1 version back
        self.assertEqual(coll['name'][0], "test_collection")
        # revert to a previous version
        self.db.revert_collection_history(id, 1)
        coll = self.db.get_collection(id)
        self.assertEqual(coll['name'][0], "test_collection")

    def test_dataset_history(self):
        # create collection
        id = self.createCollection()
        # add dataset
        self.db.add_dataset(id, "test_dataset")
        # edit dataset
        datasets = self.db.get_datasets(id)
        dataset_id = datasets[0]['uuid'][0]
        self.db.edit_dataset(dataset_id, "name", "edited_name")
        # read previous version
        dataset = self.db.read_dataset_history(dataset_id, 1) # 1 version back
        self.assertEqual(dataset['name'][0], "test_dataset")
        # revert to previous version
        self.db.revert_dataset_history(dataset_id, 1)
        dataset = self.db.get_dataset(dataset_id)
        self.assertEqual(dataset['name'][0], "test_dataset")

    def test_create_revision_on_collection(self):
        # Create a collection
        id = self.createCollection()
        # Create the revision
        rev_id = self.db.create_revision(id)
        # Get the revision
        rev = self.db.get_collection(rev_id)
        self.assertEqual(rev['name'][0], "test_collection")

    def test_edit_revision(self):
        # Create a collection
        id = self.createCollection()
        # Create a revision
        rev_id = self.db.create_revision(id)
        # Edit the revision
        self.db.edit_collection(rev_id, "name", "edited_name")
        # Get the revision
        rev = self.db.get_collection(rev_id)
        self.assertEqual(rev['name'][0], "edited_name")

    def test_publish_revision(self):
        # Create a collection
        id = self.createCollection()
        # Create a revision
        rev_id = self.db.create_revision(id)
        self.db.edit_collection(rev_id, "name", "edited_name")
        # Publish the revision
        self.db.publish_revision(rev_id)
        # Get the new collection
        coll = self.db.get_collection(id)
        self.assertEqual(coll['name'][0], "edited_name")

    def test_delete_revision(self):
        # Create a collection
        id = self.createCollection()
        # Create a revision
        rev_id = self.db.create_revision(id)
        # Delete the revision
        self.db.delete_collection(rev_id)
        # Try getting the revision
        status = self.db.get_attribute(rev_id, "published")
        self.assertEqual(status, -1)

    def test_replace_dataset_in_revision(self):
        # TODO: artifacts stuff
        pass

if __name__ == "__main__":
    unittest.main()