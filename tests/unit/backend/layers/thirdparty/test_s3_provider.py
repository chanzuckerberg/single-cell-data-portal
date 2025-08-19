import unittest
from unittest.mock import patch

from backend.layers.thirdparty.s3_exceptions import S3DeleteException
from backend.layers.thirdparty.s3_provider import S3Provider


class TestS3Provider(unittest.TestCase):
    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test__object_deletion(self, mock_client_constructor):
        test_key = "file.txt"

        with self.subTest("Raises 500 when Errors list is not empty"):
            mock_client_constructor.return_value.delete_objects.return_value = {
                "Objects": ["object"],
                "Errors": ["error"],
            }
            with self.assertRaises(S3DeleteException):
                S3Provider().delete_files("bucket", [test_key])
        with self.subTest("Does not raise errors when Errors list is empty and Objects list is not empty"):
            mock_client_constructor.return_value.delete_objects.return_value = {"Objects": ["object"], "Errors": []}
            S3Provider().delete_files("bucket", [test_key])
            self.assertTrue(True)  # No exception raised

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_delete_files_empty_keys(self, mock_client_constructor):
        """Test that delete_files handles empty key lists correctly"""
        mock_client_constructor.return_value.delete_objects.return_value = {"Objects": [], "Errors": []}
        s3_provider = S3Provider()

        # Test with empty list - should not call delete_objects
        s3_provider.delete_files("bucket", [])
        mock_client_constructor.return_value.delete_objects.assert_not_called()

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_delete_files_with_keys(self, mock_client_constructor):
        """Test that delete_files works with valid keys"""
        mock_client_constructor.return_value.delete_objects.return_value = {"Objects": ["object"], "Errors": []}

        s3_provider = S3Provider()
        keys = ["file1.txt", "file2.txt"]

        s3_provider.delete_files("bucket", keys)

        # Should call delete_objects with all keys
        mock_client_constructor.return_value.delete_objects.assert_called_once_with(
            Bucket="bucket", Delete={"Objects": [{"Key": "file1.txt"}, {"Key": "file2.txt"}]}
        )

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_delete_files_batch_slicing_fix(self, mock_client_constructor):
        """Test that batch slicing works correctly with the fix"""
        mock_client_constructor.return_value.delete_objects.return_value = {"Objects": ["object"], "Errors": []}

        # Create more keys than batch size to test slicing
        large_key_list = [f"key_{i}.txt" for i in range(1500)]  # More than AWS_S3_MAX_ITEMS_PER_BATCH (1000)

        s3_provider = S3Provider()
        s3_provider.delete_files("bucket", large_key_list)

        # Should be called twice: once for first 1000, once for remaining 500
        self.assertEqual(mock_client_constructor.return_value.delete_objects.call_count, 2)


if __name__ == "__main__":
    unittest.main()
