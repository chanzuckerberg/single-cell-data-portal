import unittest
from unittest.mock import patch

from backend.layers.thirdparty.s3_exceptions import S3DeleteException, IllegalS3RecursiveDelete
from backend.layers.thirdparty.s3_provider import S3Provider


class TestS3Provider(unittest.TestCase):
    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test__object_deletion(self, mock_client_constructor):
        safe_key = "12345678-1234-5678-9012-123456789012/file.txt"
        
        with self.subTest("Raises 500 when Errors list is not empty"):
            mock_client_constructor.return_value.delete_objects.return_value = {
                "Objects": ["object"],
                "Errors": ["error"],
            }
            with self.assertRaises(S3DeleteException):
                S3Provider().delete_files("bucket", [safe_key])
        with self.subTest("Does not raise errors when Errors list is empty and Objects list is not empty"):
            mock_client_constructor.return_value.delete_objects.return_value = {"Objects": ["object"], "Errors": []}
            S3Provider().delete_files("bucket", [safe_key])
            self.assertTrue(True)  # No exception raised

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_delete_files_empty_keys(self, mock_client_constructor):
        """Test that delete_files handles empty key lists correctly"""
        s3_provider = S3Provider()
        
        # Test with empty list
        s3_provider.delete_files("bucket", [])
        mock_client_constructor.return_value.delete_objects.assert_not_called()
        
        # Test with list of empty strings
        s3_provider.delete_files("bucket", ["", "  ", None])
        mock_client_constructor.return_value.delete_objects.assert_not_called()

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_delete_files_filters_invalid_keys(self, mock_client_constructor):
        """Test that delete_files filters out empty/invalid keys"""
        mock_client_constructor.return_value.delete_objects.return_value = {
            "Objects": ["object"], "Errors": []
        }
        
        s3_provider = S3Provider()
        valid_key1 = "12345678-1234-5678-9012-123456789012/file1.txt"
        valid_key2 = "87654321-4321-8765-2109-876543210987/file2.txt"
        
        s3_provider.delete_files("bucket", [valid_key1, "", "  ", valid_key2])
        
        # Should only call delete_objects with valid keys
        mock_client_constructor.return_value.delete_objects.assert_called_once_with(
            Bucket="bucket",
            Delete={"Objects": [{"Key": valid_key1}, {"Key": valid_key2}]}
        )

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_delete_files_batch_slicing_fix(self, mock_client_constructor):
        """Test that batch slicing works correctly with the fix"""
        mock_client_constructor.return_value.delete_objects.return_value = {
            "Objects": ["object"], "Errors": []
        }
        
        # Create more keys than batch size to test slicing - use properly prefixed keys
        large_key_list = [f"12345678-1234-5678-9012-123456789012/key_{i}" for i in range(1500)]  # More than AWS_S3_MAX_ITEMS_PER_BATCH (1000)
        
        s3_provider = S3Provider()
        s3_provider.delete_files("bucket", large_key_list)
        
        # Should be called twice: once for first 1000, once for remaining 500
        self.assertEqual(mock_client_constructor.return_value.delete_objects.call_count, 2)

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_delete_files_blocks_dangerous_keys(self, mock_client_constructor):
        """Test that delete_files blocks deletion of root-level or insufficiently prefixed objects"""
        s3_provider = S3Provider()
        
        # Test root-level files (no slash)
        with self.assertRaises(IllegalS3RecursiveDelete):
            s3_provider.delete_files("bucket", ["file1.txt", "file2.txt"])
        
        # Test insufficiently prefixed files (prefix < 8 chars)
        with self.assertRaises(IllegalS3RecursiveDelete):
            s3_provider.delete_files("bucket", ["abc/file.txt", "short/file.txt"])
        
        # Test mixed dangerous and safe keys - should still raise exception
        with self.assertRaises(IllegalS3RecursiveDelete):
            s3_provider.delete_files("bucket", ["12345678-1234/safe.txt", "unsafe.txt"])
        
        # Verify delete_objects was never called for dangerous operations
        mock_client_constructor.return_value.delete_objects.assert_not_called()

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_delete_files_allows_safe_keys(self, mock_client_constructor):
        """Test that delete_files allows properly prefixed objects"""
        mock_client_constructor.return_value.delete_objects.return_value = {
            "Objects": ["object"], "Errors": []
        }
        
        s3_provider = S3Provider()
        
        # These should be allowed (8+ character prefix before slash)
        safe_keys = [
            "12345678-1234-5678-9012-123456789012/file1.txt",
            "longenoughprefix/subfolder/file2.txt",
            "uuid-like-prefix/data/file3.txt"
        ]
        
        # Should not raise exception
        s3_provider.delete_files("bucket", safe_keys)
        
        # Verify delete_objects was called
        mock_client_constructor.return_value.delete_objects.assert_called_once()


if __name__ == "__main__":
    unittest.main()
