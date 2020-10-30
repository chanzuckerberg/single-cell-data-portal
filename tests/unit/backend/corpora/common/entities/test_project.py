import unittest
from datetime import datetime

from sqlalchemy.exc import SQLAlchemyError

from backend.corpora.common.corpora_orm import ProjectLinkType, DbProjectLink, CollectionVisibility, DbDataset
from backend.corpora.common.entities import Dataset
from backend.corpora.common.entities.project import Project
from backend.corpora.common.utils.db_utils import DbUtils
from tests.unit.backend.utils import BogusProjectParams, BogusDatasetParams


class TestProject(unittest.TestCase):
    def setUp(self):
        self.uuid = "test_project_id"
        self.visibility = CollectionVisibility.PUBLIC.name
        self.db = DbUtils()

    def tearDown(self):
        self.db.session.rollback()
        self.db.close()

    def test__get__ok(self):
        key = (self.uuid, self.visibility)

        project = Project.get(key)

        # Verify Columns
        self.assertEqual(project.name, "test_project")
        self.assertEqual(project.owner, "test_user_id")

        # Verify Dataset relationship
        dataset = project.datasets[0]
        self.assertIsInstance(dataset, DbDataset)
        self.assertEqual(dataset.id, "test_dataset_id")
        self.assertEqual(len(dataset.assay), 1)
        self.assertDictEqual(dataset.assay[0], {"ontology_term_id": "test_obo", "label": "test_assay"})

        # Verify Link relationship
        self.assertIsInstance(project.links[0], DbProjectLink)
        self.assertEqual(project.links[0].id, "test_project_link_id")

    def test__get__does_not_exist(self):
        non_existent_key = ("non_existent_id", self.visibility)

        self.assertEqual(Project.get(non_existent_key), None)

    def test__get__invalid_visibility(self):
        invalid_visibility_key = (self.uuid, "invalid_visibility")
        with self.assertRaises(SQLAlchemyError):
            Project.get(invalid_visibility_key)

    def test__create__ok(self):
        """
        Create a project with a variable number of links.
        """

        link_params = {"link_url": "fake_url", "link_type": ProjectLinkType.PROTOCOL.name}
        project_params = BogusProjectParams.get()

        for i in range(3):
            with self.subTest(i):
                project = Project.create(links=[link_params] * i, **project_params)

                project_key = (project.id, project.visibility)
                expected_links = project.links

                # Expire all local object and retrieve them from the DB to make sure the transactions went through.
                Project.db.session.expire_all()

                actual_project = Project.get(project_key)
                self.assertEqual(project_key, (actual_project.id, actual_project.visibility))
                self.assertCountEqual(expected_links, actual_project.links)

    def test__list_in_time_range__ok(self):
        created_before = Project.create(**BogusProjectParams.get(), created_at=datetime.fromtimestamp(10))
        from_date = 20
        created_inbetween = Project.create(**BogusProjectParams.get(), created_at=datetime.fromtimestamp(30))
        to_date = 40
        created_after = Project.create(**BogusProjectParams.get(), created_at=datetime.fromtimestamp(50))

        with self.subTest("from_date"):
            # Projects from_date are returned.
            actual_projects = Project.list_attributes_in_time_range(from_date=from_date)
            self.assertTrue(all([p["created_at"].timestamp() > from_date for p in actual_projects]))
            expected_ids = [created_inbetween.id, created_after.id, "test_project_id"]
            actual_ids = [p["id"] for p in actual_projects]
            # Check if the test ids we created are present.
            # As a result of other tests, more projects have likely been created and will be return in the results,
            # so we can't do an exact match.
            self.assertTrue(set(expected_ids).issubset(actual_ids))

        with self.subTest("to_date"):
            # Projects to_date are returned.
            actual_projects = Project.list_attributes_in_time_range(to_date=to_date)
            self.assertTrue(all([p["created_at"].timestamp() < to_date for p in actual_projects]))
            expected_ids = [created_before.id, created_inbetween.id]
            actual_ids = [p["id"] for p in actual_projects]
            self.assertCountEqual(expected_ids, actual_ids)

        with self.subTest("from_date->to_date"):
            # Projects between to_date and from_date are returned.
            actual_projects = Project.list_attributes_in_time_range(to_date=to_date, from_date=from_date)
            self.assertTrue(all([p["created_at"].timestamp() > from_date for p in actual_projects]))
            self.assertTrue(all([p["created_at"].timestamp() < to_date for p in actual_projects]))
            expected_ids = [created_inbetween.id]
            actual_ids = [p["id"] for p in actual_projects]
            self.assertCountEqual(expected_ids, actual_ids)

        with self.subTest("No parameters"):
            """All projects are returned."""
            actual_projects = Project.list_attributes_in_time_range()
            expected_ids = [created_before.id, created_inbetween.id, created_after.id]
            actual_ids = [p["id"] for p in actual_projects]
            # Check if the test ids we created are present.
            # As a result of other tests, more projects have likely been created and will be return in the results,
            # so we can't do an exact match.
            self.assertTrue(set(expected_ids).issubset(actual_ids))

    def test__list_projects_in_time_range___ok(self):
        """Public projects are returned"""
        from_date = 10
        created_at = datetime.fromtimestamp(20)
        to_date = 30
        expected_project = Project.create(
            **BogusProjectParams.get(created_at=created_at, visibility=CollectionVisibility.PUBLIC.name)
        )

        # Generate a project with Private visibility. This should not show up in the results.
        Project.create(**BogusProjectParams.get(created_at=created_at, visibility=CollectionVisibility.PRIVATE.name))

        actual_projects = Project.list_projects_in_time_range(to_date=to_date, from_date=from_date)
        expected_projects = [{"created_at": created_at, "id": expected_project.id}]
        self.assertCountEqual(expected_projects, actual_projects)

    def test__list__ok(self):
        generate = 2
        generated_ids = [Project.create(**BogusProjectParams.get()).id for _ in range(generate)]
        projects = Project.list()
        self.assertTrue(set(generated_ids).issubset([p.id for p in projects]))

    def test__cascade_delete_project__ok(self):
        test_project = Project.create(**BogusProjectParams.get(links=[{}]))
        db = DbUtils()
        test_project_id = test_project.id
        test_link_id = test_project.links[0].id

        # Check if the link exists.
        expected_id = test_link_id
        actual_id = db.query([DbProjectLink], [DbProjectLink.id == test_link_id])[0].id
        self.assertEqual(expected_id, actual_id)

        # The Project is deleted
        db.session.delete(test_project.db_object)
        expected_results = None
        actual_results = Project.get_project(test_project_id)
        self.assertEqual(expected_results, actual_results)

        # The link should also be deleted.
        expected_results = []
        actual_results = db.query([DbProjectLink], [DbProjectLink.id == test_link_id])
        self.assertEqual(expected_results, actual_results)

    def test__cascade_delete_project_with_dataset__ok(self):
        db = DbUtils()
        test_project = Project.create(**BogusProjectParams.get())
        expected_project_id = test_project.id
        test_dataset = Dataset.create(
            **BogusDatasetParams.get(collection_id=test_project.id, collection_visibility=test_project.visibility)
        )
        expected_dataset_id = test_dataset.id

        # The Project is deleted
        db.session.delete(test_project.db_object)
        db.session.expire_all()
        expected_results = None
        actual_results = Project.get_project(expected_project_id)
        self.assertEqual(expected_results, actual_results)

        # The dataset is delete
        expected_results = None
        actual_results = Dataset.get(expected_dataset_id)
        self.assertEqual(expected_results, actual_results)
