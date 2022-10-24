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
import copy
from unittest.mock import Mock
from backend.corpora.common.providers.crossref_provider import CrossrefException, CrossrefProvider
from backend.layers.business.business import BusinessLogic, UserInfo
from backend.layers.common.entities import CollectionMetadata, Link
from backend.layers.thirdparty.crossref_provider import CrossrefProviderInterface
from tests.unit.backend.layers.persistence.persistence_mock import DatabaseProviderMock

# test fixtures - TODO move them to a fixture

metadata1 = CollectionMetadata(
    "test collection 1", 
    "description of test collection 1",
    "test_user_1",
    "scientist",
    "scientist@czi.com",
    []
)

user_info = UserInfo(
    "test_user_1",
    "fake token"
)

# end fixtures

class BaseBusinessLogicTestCase(unittest.TestCase):

    def setUp(self) -> None:
        self.database_provider = DatabaseProviderMock()

        # By default does nothing. Can be mocked by single test cases.
        self.crossref_provider = CrossrefProviderInterface()

        self.business_logic = BusinessLogic(
            database_provider=self.database_provider, 
            crossref_provider=self.crossref_provider
        )
        return super().setUp()


class TestCreateCollection(BaseBusinessLogicTestCase):

    def create_collection_ok(self):
        cv = self.business_logic.create_collection(metadata1, user_info)
        cv_from_database = self.database_provider.get_collection_version(cv.version_id)
        self.assertEqual(cv, cv_from_database)

    def create_collection_with_links_ok(self):
        good_links = [
            Link("test link 1", "protocol", "http://example.com/protocol"),
            Link("test link 2", "other", "http://example.com/other"),
        ]
        metadata = copy.deepcopy(metadata1)
        metadata.links = good_links
        cv = self.business_logic.create_collection(metadata, user_info)
        cv_from_database = self.database_provider.get_collection_version(cv.version_id)
        self.assertEqual(cv.metadata.links, cv_from_database.metadata.links)

    def create_collection_with_bad_links_fail(self):
        bad_links = [
            Link("test bad link", "other", "incorrect_url")
        ]
        metadata = copy.deepcopy(metadata1)
        metadata.links = bad_links
        # TODO: Create InvalidLinkException
        # TODO: the current method has error message validation - implement a more sophisticate version
        with self.assertRaises(InvalidLinkException):
            self.business_logic.create_collection(metadata, user_info)

    def create_collection_with_valid_doi_ok(self):
        links_with_doi = [
            Link("test doi", "doi", "http://good.doi")
        ]
        metadata = copy.deepcopy(metadata1)
        metadata.links = links_with_doi

        expected_publiser_metadata = {"authors": ["Test Author"]}

        self.crossref_provider.fetch_metadata = Mock(return_value=expected_publiser_metadata)

        cv = self.business_logic.create_collection(metadata, user_info)

        self.crossref_provider.fetch_metadata.assert_called_with("http://good.doi")

        cv_from_database = self.database_provider.get_collection_version(cv.version_id)
        self.assertEqual(1, len(cv_from_database.metadata.links))
        self.assertEqual(cv_from_database.metadata.links[0].uri, "http://good.doi")
        self.assertIsNotNone(cv_from_database.publisher_metadata)
        self.assertEqual(cv_from_database.publisher_metadata, expected_publiser_metadata)


    def create_collection_with_invalid_doi_fail(self):
        links_with_doi = [
            Link("test doi", "doi", "http://bad.doi")
        ]

        metadata = copy.deepcopy(metadata1)
        metadata.links = links_with_doi

        # TODO: make sure that we don't need different actions depending on which exception
        self.crossref_provider.fetch_metadata = Mock(side_effect=CrossrefException("Error!"))

        cv = self.business_logic.create_collection(metadata, user_info)
        # TODO: create CollectionCreationException
        with self.assertRaises(CollectionCreationException):
            self.business_logic.create_collection(metadata, user_info)

    def create_collection_unauthorized_fail(self):
        # TODO: AUTHORIZATION
        pass

class TestDeleteCollection(BaseBusinessLogicTestCase):

    def delete_collection_ok(self):
        # TODO: do we even need collection deletion? I don't think this is supported.
        pass

    def delete_collection_unauthorized_fail(self):
        # TODO: AUTHORIZATION
        pass

class TestGetCollections(BaseBusinessLogicTestCase):

    def get_collection_ok(self):
        self.database_provider.create_collection
        created_version = self.business_logic.create_collection(metadata1, user_info)
        self.business_logic.get_collection()


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

class TestUpdateCollection(BaseBusinessLogicTestCase):

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

class TestCreateCollectionVersion(BaseBusinessLogicTestCase):

    def create_collection_version_ok(self):
        pass

    def create_collection_version_unauthorized_fail(self):
        pass

class TestDeleteCollectionVersion(BaseBusinessLogicTestCase):

    def delete_collection_version_ok(self):
        pass

    def delete_collection_version_unauthorized_fail(self):
        pass

class TestAddDataset(BaseBusinessLogicTestCase):

    def add_dataset_to_collection_ok(self):
        pass

    def add_dataset_to_collection_fail(self):
        # TODO: are there failures to be tested here?
        pass

class TestUpdateRemoveDataset(BaseBusinessLogicTestCase):

    def remove_dataset_from_collection_ok(self):
        pass

    def replace_dataset_in_collection_ok(self):
        pass

class TestPublishCollectionVersion(BaseBusinessLogicTestCase):

    def publish_version_ok(self):
        pass

    def publish_version_unauthorized_fail(self):
        pass

class TestGetDataset(BaseBusinessLogicTestCase):

    def get_all_datasets_ok(self):
        pass

    def get_dataset_artifacts_ok(self):
        pass

class TestDatasetStatus(BaseBusinessLogicTestCase):

    def update_dataset_status_ok(self):
        pass

    def get_dataset_status_ok(self):
        pass


if __name__ == '__main__':
    unittest.main()

