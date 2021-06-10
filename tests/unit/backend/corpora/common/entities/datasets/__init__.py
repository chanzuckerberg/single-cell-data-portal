import typing

from backend.corpora.common.corpora_orm import (
    Base,
    DatasetArtifactFileType,
    DatasetArtifactType,
    UploadStatus,
    ValidationStatus,
)
from tests.unit.backend.fixtures.mock_aws_test_case import CorporaTestCaseUsingMockAWS
from tests.unit.backend.utils import BogusDatasetParams


class TestDataset(CorporaTestCaseUsingMockAWS):
    def setUp(self):
        super().setUp()
        self.uuid = "test_dataset_id"
        self.bucket_name = self.CORPORA_TEST_CONFIG["bucket_name"]

    def create_dataset_with_artifacts(self, artifact_count=1, artifact_params=None):
        """
        Create a dataset with a variable number of artifacts
        """
        if not artifact_params:
            artifact_params = dict(
                filename="filename_x",
                filetype=DatasetArtifactFileType.H5AD,
                type=DatasetArtifactType.ORIGINAL,
                user_submitted=True,
                s3_uri="some_uri",
            )

        dataset_params = BogusDatasetParams.get()
        dataset = self.generate_dataset(
            self.session,
            **dataset_params,
            artifacts=[artifact_params] * artifact_count,
            processing_status={
                "upload_progress": 9 / 13,
                "upload_status": UploadStatus.UPLOADING,
                "validation_status": ValidationStatus.NA,
            },
        )
        return dataset

    def assertRowsDeleted(self, tests: typing.List[typing.Tuple[str, Base]]):
        """
        Verify if rows have been deleted from the database.
        :param tests: a list of tuples with (primary_key, table)
        """
        self.session.expire_all()
        for p_key, table in tests:
            if len(p_key) == 2:
                # handle the special case for collections with a composite primary key
                actual = [
                    i for i in self.session.query(table).filter(table.id == p_key[0], table.visibility == p_key[1])
                ]
            else:
                actual = [i for i in self.session.query(table).filter(table.id == p_key)]
            self.assertFalse(actual, f"Row not deleted {table.__name__}:{p_key}")

    def assertRowsExist(self, tests):
        """
        Verify if rows exist in the database.
        :param tests: a list of tuples with (primary_key, table)
        """
        self.session.expire_all()
        for p_key, table in tests:
            if len(p_key) == 2:
                # handle the special case for collections with a composite primary key
                actual = [
                    i for i in self.session.query(table).filter(table.id == p_key[0], table.visibility == p_key[1])
                ]
            else:
                actual = [i for i in self.session.query(table).filter(table.id == p_key)]
            self.assertTrue(actual, f"Row does not exist {table.__name__}:{p_key}")
