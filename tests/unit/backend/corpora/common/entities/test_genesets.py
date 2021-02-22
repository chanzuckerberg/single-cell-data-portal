from tests.unit.backend.fixtures.data_portal_test_case import DataPortalTestCase


class TestGeneSets(DataPortalTestCase):
    def setUp(self):
        super().setUp()

    def test_create_and_retrieve_a_geneset(self):
        test_uuid = "test_genest_id"
        geneset = self.generate_geneset(self.session, id=test_uuid)
        self.assertEqual(test_uuid, geneset.id)

    def test_retrieve_a_geneset_that_does_not_exist(self):
        pass

    def test_retrieve_list_of_genesets_for_colledtion(self):
        pass

    def test_link_geneset_to_dataset(self):
        pass

    def test_retrieve_all_genesets_for_a_dataset(self):
        pass

    def test_retrieve_all_datasets_for_a_geneset(self):
        pass
