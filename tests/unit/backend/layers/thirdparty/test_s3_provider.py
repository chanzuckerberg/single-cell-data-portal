import unittest
from unittest.mock import patch

from backend.common.utils.http_exceptions import ServerErrorHTTPException
from backend.layers.thirdparty.s3_provider import S3Provider


class TestS3Provider(unittest.TestCase):
    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test__object_deletion(self, mock_client_constructor):
        with self.subTest("Raises 500 when Errors list is not empty"):
            mock_client_constructor.return_value.delete_objects.return_value = {
                "Objects": ["object"],
                "Errors": ["error"],
            }
            with self.assertRaises(ServerErrorHTTPException):
                S3Provider().delete_files("bucket", ["key"])
        with self.subTest("Does not raise errors when Errors list is empty and Objects list is not empty"):
            mock_client_constructor.return_value.delete_objects.return_value = {"Objects": ["object"], "Errors": []}
            S3Provider().delete_files("bucket", ["key"])
            self.assertTrue(True)  # No exception raised


if __name__ == "__main__":
    unittest.main()
