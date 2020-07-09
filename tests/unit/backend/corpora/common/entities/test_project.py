import logging
import time
import unittest
from datetime import datetime

from backend.corpora.common.corpora_orm import (
    ProjectLinkType,
    DbProjectLink,
    ProcessingState,
    ValidationState,
    ProjectStatus,
    DbDataset,
    DbUser,
)
from backend.corpora.common.entities.entity import logger as entity_logger
from backend.corpora.common.entities.project import Project
from tests.unit.backend.corpora.common.entities import get_ids


class ProjectParams:
    i = 0

    @classmethod
    def get(cls):
        cls.i += 1
        return dict(
            status=ProjectStatus.EDIT.name,
            name=f"Created Project {cls.i}",
            description="test",
            owner="test_user_id",
            s3_bucket="s3://fakebucket",
            tc_uri="https://fakeurl",
            needs_attestation=False,
            processing_state=ProcessingState.IN_VALIDATION.name,
            validation_state=ValidationState.NOT_VALIDATED.name,
        )


class TestProject(unittest.TestCase):
    def setUp(self):
        self.uuid = "test_project_id"
        self.status = ProjectStatus.LIVE.name

    def test__get__ok(self):
        key = (self.uuid, self.status)

        project = Project.get(key)

        # Verify Columns
        self.assertEqual(project.name, "test_project")
        self.assertEqual(project.owner, "test_user_id")

        # Verify User relationship
        self.assertIsInstance(project.user, DbUser)
        self.assertEqual(project.user.id, "test_user_id")

        # Verify Dataset relationship
        dataset = project.datasets[0]
        self.assertIsInstance(dataset, DbDataset)
        self.assertEqual(dataset.id, "test_dataset_id")
        self.assertEqual(dataset.assay, "test_assay")

        # Verify Link relationship
        self.assertIsInstance(project.links[0], DbProjectLink)
        self.assertEqual(project.links[0].id, "test_project_link_id")

    def test__get__does_not_exist(self):
        non_existent_key = ("non_existent_id", self.status)

        self.assertEqual(Project.get(non_existent_key), None)

    def test__get__invalid_status(self):
        invalid_status_key = (self.uuid, "invalid_status")
        with self.assertLogs(entity_logger, logging.INFO) as logs:
            self.assertEqual(Project.get(invalid_status_key), None)
        self.assertIn("Unable to find a row with primary key", logs.output[0])
        self.assertEqual(Project.get(invalid_status_key), None)

    def test__create__ok(self):
        """
        Create a project with a variable number of links.
        """

        link_params = {"link_url": "fake_url", "link_type": ProjectLinkType.PROTOCOL.name}
        project_params = ProjectParams.get()

        for i in range(3):
            with self.subTest(i):
                project = Project.create(**project_params, links=[link_params] * i)

                project_key = (project.id, project.status)
                link_ids = get_ids(project.links)

                # Expire all local object and retireve them from the DB to make sure the transactions went through.
                Project.db.session.expire_all()

                project = Project.get(project_key)
                self.assertIsNotNone(project)
                self.assertEqual(project_key, (project.id, project.status))
                self.assertEqual(link_ids, get_ids(project.links))

        with self.subTest("With existing rows"):
            project = Project.create(**project_params, links=[{"id": "test_project_link_id"}])

            project_key = (project.id, project.status)
            link_ids = get_ids(project.links)

            # Expire all local object and retrieve them from the DB to make sure the transactions went through.
            Project.db.session.expire_all()

            project = Project.get(project_key)
            self.assertIsNotNone(project)
            self.assertEqual(project_key, (project.id, project.status))

            self.assertEqual(link_ids, get_ids(project.links))
            self.assertNotEqual(["test_project_link_id"], get_ids(project.links))

    def test__list_in_time_range__ok(self):
        generate = 5
        sleep = 0.5
        from_ids = []
        to_ids = []

        from_date = datetime.now().timestamp()
        time.sleep(sleep)
        for i in range(generate):
            from_ids.append(Project.create(**ProjectParams.get()).id)

        with self.subTest("Test from_date"):
            # Projects from_date are returned.
            projects = Project.list_in_time_range(from_date=from_date)
            for project in projects:
                self.assertGreaterEqual(project["created_at"].timestamp(), from_date)

        to_date = datetime.now().timestamp()
        time.sleep(sleep)
        for i in range(generate):
            to_ids.append(Project.create(**ProjectParams.get()).id)

        with self.subTest("Test to_date and from_date"):
            # Projects between to_date and from_date are returned.
            projects = Project.list_in_time_range(to_date=to_date, from_date=from_date)
            for project in projects:
                self.assertGreaterEqual(project["created_at"].timestamp(), from_date)
                self.assertLessEqual(project["created_at"].timestamp(), to_date)

        with self.subTest("Test to_date"):
            # All created projects are returned.
            projects = Project.list_in_time_range(to_date=to_date)
            self.assertGreater(len(projects), generate)
            for project in projects:
                self.assertLessEqual(project["created_at"].timestamp(), to_date)

        with self.subTest("Test no date"):
            projects = Project.list_in_time_range()
            test_ids = set([*to_ids, *from_ids])
            project_ids = set([p["id"] for p in projects])
            self.assertTrue(test_ids.issubset(project_ids))

        with self.subTest("Verify columns"):
            self.assertIsInstance(projects[0]["id"], str)
            self.assertIsInstance(projects[0]["created_at"], datetime)
            self.assertIsInstance(projects[0]["name"], str)

    def test__list__ok(self):
        generate = 5

        for i in range(generate):
            Project.create(**ProjectParams.get())

        projects = Project.list()
        self.assertGreaterEqual(len(projects), generate)
        self.assertTrue(all([isinstance(i, Project) for i in projects]))
