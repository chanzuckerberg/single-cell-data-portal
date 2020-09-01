from tests.unit.backend.corpora import CorporaTestCaseUsingMockAWS
from backend.corpora.common.utils.s3_utils import generate_file_url, head_file


class TestS3Utils(CorporaTestCaseUsingMockAWS):
    def test_generate_file_url(self):
        file_prefix = "test_prefix"
        bucket_name = self.CORPORA_TEST_CONFIG["bucket_name"]

        with self.subTest("Generate URL OK"):
            url = generate_file_url(bucket_name, file_prefix)
            self.assertIn(file_prefix, url)
            self.assertIn("Expires=", url)

    def test_head_file(self):
        file_name = "fake_file.txt"
        content = "hello world!"
        bucket = self.CORPORA_TEST_CONFIG["bucket_name"]
        self.create_s3_object(file_name, bucket, content=content)

        response = head_file(bucket, file_name)
        self.assertEqual(len(content), response["ContentLength"])

    def test_head_file__not_found(self):
        file_name = "fake_file.txt"
        bucket_name = self.CORPORA_TEST_CONFIG["bucket_name"]
        self.assertIs(None, head_file(bucket_name, file_name))
