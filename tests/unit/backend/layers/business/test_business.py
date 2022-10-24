# from datetime import datetime

# from backend.corpora.common.corpora_orm import CollectionLinkType, DbCollectionLink, CollectionVisibility, DbDataset
# from backend.corpora.common.entities import Dataset
# from backend.corpora.common.entities.collection import Collection
# from backend.corpora.common.entities.geneset import Geneset, GenesetDatasetLink
# from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase
# from tests.unit.backend.utils import BogusCollectionParams, BogusDatasetParams



    # def setUp(self) -> None:
    #     self.database_provider = DatabaseProviderMock()
    #     return super().setUp()



import unittest
from tests.unit.backend.layers.persistence.persistence_mock import DatabaseProviderMock

class TestCreateCollection(unittest.TestCase):

    def create_collection_ok(self):
        pass

    def create_collection_with_links_ok(self):
        pass

    def create_collection_with_bad_links_fail(self):
        pass

    def create_collection_with_valid_doi_ok(self):
        pass

    def create_collection_with_invalid_doi_fail(self):
        pass

    def create_collection_unauthorized_fail(self):
        pass

class TestDeleteCollection(unittest.TestCase):

    def delete_collection_ok(self):
        pass

    def delete_collection_unauthorized_fail(self):
        pass

class TestGetCollections(unittest.TestCase):

    def get_collection_ok(self):
        pass

    def get_collection_missing_fail(self):
        pass

    def get_all_collections_ok(self):
        pass

    def get_all_collections_with_owner_ok(self):
        pass

    def get_all_collections_with_date_range_ok(self):
        pass

    def get_all_collections_visibility_ok(self):
        pass

class TestUpdateCollection(unittest.TestCase):

    def update_collection_ok(self):
        pass

    def update_collection_same_doi(self):
        pass

    def update_collection_change_doi(self):
        pass

    def update_collection_wrong_params_fail(self):
        pass

    def update_collection_unauthorized_fail(self):
        pass

class TestCreateCollectionVersion(unittest.TestCase):

    def create_collection_version_ok(self):
        pass

    def create_collection_version_unauthorized_fail(self):
        pass

class TestDeleteCollectionVersion(unittest.TestCase):

    def delete_collection_version_ok(self):
        pass

    def delete_collection_version_unauthorized_fail(self):
        pass

class TestAddDataset(unittest.TestCase):

    def add_dataset_to_collection_ok(self):
        pass

    def add_dataset_to_collection_fail(self):
        # TODO: are there failures to be tested here?
        pass

class TestUpdateRemoveDataset(unittest.TestCase):

    def remove_dataset_from_collection_ok(self):
        pass

    def replace_dataset_in_collection_ok(self):
        pass

class TestPublishCollectionVersion(unittest.TestCase):

    def publish_version_ok(self):
        pass

    def publish_version_unauthorized_fail(self):
        pass

class TestGetDataset(unittest.TestCase):

    def get_all_datasets_ok(self):
        pass

    def get_dataset_artifacts_ok(self):
        pass

class TestDatasetStatus(unittest.TestCase):

    def update_dataset_status_ok(self):
        pass

    def get_dataset_status_ok(self):
        pass


if __name__ == '__main__':
    unittest.main()

