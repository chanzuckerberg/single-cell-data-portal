import os
import tempfile
import unittest
from unittest.mock import MagicMock, patch

from backend.layers.thirdparty.s3_exceptions import S3DeleteException
from backend.layers.thirdparty.s3_provider import SIZE_50GB, S3Provider


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


class TestS3ProviderCacheControl(unittest.TestCase):
    """Test Cache-Control metadata setting for large files (> 50GB)"""

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_upload_small_file_no_cache_control(self, mock_client_constructor):
        """Small files (< 50GB) should NOT have Cache-Control metadata added"""
        mock_client = MagicMock()
        mock_client_constructor.return_value = mock_client

        s3_provider = S3Provider()

        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            # Create small file (1 MB)
            tmp_file.write(b"x" * (1024 * 1024))
            tmp_file.flush()
            file_path = tmp_file.name

        try:
            extra_args = {"ACL": "bucket-owner-full-control"}

            s3_provider.upload_file(file_path, "test-bucket", "test-key", extra_args)

            # Verify Cache-Control was NOT added
            call_args = mock_client.upload_file.call_args
            self.assertNotIn("CacheControl", call_args[1]["ExtraArgs"])
            self.assertEqual(call_args[1]["ExtraArgs"]["ACL"], "bucket-owner-full-control")

        finally:
            os.unlink(file_path)

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    @patch("os.path.getsize")
    def test_upload_large_file_sets_cache_control(self, mock_getsize, mock_client_constructor):
        """Large files (> 50GB) should have Cache-Control set to no-cache"""
        mock_client = MagicMock()
        mock_client_constructor.return_value = mock_client

        s3_provider = S3Provider()

        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            tmp_file.write(b"test data")
            tmp_file.flush()
            file_path = tmp_file.name

        try:
            # Mock file size to be > 50GB
            mock_getsize.return_value = SIZE_50GB + 1

            extra_args = {"ACL": "bucket-owner-full-control"}

            s3_provider.upload_file(file_path, "test-bucket", "test-key", extra_args)

            # Verify Cache-Control was added
            call_args = mock_client.upload_file.call_args
            self.assertEqual(call_args[1]["ExtraArgs"]["CacheControl"], "no-cache, no-store, must-revalidate")
            self.assertEqual(call_args[1]["ExtraArgs"]["ACL"], "bucket-owner-full-control")

        finally:
            os.unlink(file_path)

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    @patch("os.path.getsize")
    def test_upload_exactly_50gb_no_cache_control(self, mock_getsize, mock_client_constructor):
        """File exactly at 50GB threshold should NOT have Cache-Control (boundary test)"""
        mock_client = MagicMock()
        mock_client_constructor.return_value = mock_client

        s3_provider = S3Provider()

        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            tmp_file.write(b"test data")
            tmp_file.flush()
            file_path = tmp_file.name

        try:
            # Mock file size to be exactly 50GB
            mock_getsize.return_value = SIZE_50GB

            extra_args = {"ACL": "bucket-owner-full-control"}

            s3_provider.upload_file(file_path, "test-bucket", "test-key", extra_args)

            # Verify Cache-Control was NOT added (only > 50GB, not >=)
            call_args = mock_client.upload_file.call_args
            self.assertNotIn("CacheControl", call_args[1]["ExtraArgs"])

        finally:
            os.unlink(file_path)

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    @patch("os.path.getsize")
    def test_upload_large_file_preserves_extra_args(self, mock_getsize, mock_client_constructor):
        """Cache-Control should be added without removing other ExtraArgs"""
        mock_client = MagicMock()
        mock_client_constructor.return_value = mock_client

        s3_provider = S3Provider()

        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            tmp_file.write(b"test data")
            tmp_file.flush()
            file_path = tmp_file.name

        try:
            mock_getsize.return_value = SIZE_50GB + 1

            extra_args = {
                "ACL": "bucket-owner-full-control",
                "ContentType": "application/octet-stream",
                "Metadata": {"custom": "value"},
            }

            s3_provider.upload_file(file_path, "test-bucket", "test-key", extra_args)

            call_args = mock_client.upload_file.call_args
            uploaded_args = call_args[1]["ExtraArgs"]

            # Verify all original args preserved
            self.assertEqual(uploaded_args["ACL"], "bucket-owner-full-control")
            self.assertEqual(uploaded_args["ContentType"], "application/octet-stream")
            self.assertEqual(uploaded_args["Metadata"], {"custom": "value"})
            # And Cache-Control added
            self.assertEqual(uploaded_args["CacheControl"], "no-cache, no-store, must-revalidate")

        finally:
            os.unlink(file_path)

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    @patch("os.path.getsize")
    def test_upload_large_file_with_none_extra_args(self, mock_getsize, mock_client_constructor):
        """Should handle None extra_args gracefully for large files"""
        mock_client = MagicMock()
        mock_client_constructor.return_value = mock_client

        s3_provider = S3Provider()

        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            tmp_file.write(b"test data")
            tmp_file.flush()
            file_path = tmp_file.name

        try:
            mock_getsize.return_value = SIZE_50GB + 1

            s3_provider.upload_file(file_path, "test-bucket", "test-key", None)

            call_args = mock_client.upload_file.call_args
            uploaded_args = call_args[1]["ExtraArgs"]

            # Should only have Cache-Control
            self.assertEqual(uploaded_args["CacheControl"], "no-cache, no-store, must-revalidate")

        finally:
            os.unlink(file_path)

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    @patch("os.path.getsize")
    def test_upload_large_file_does_not_modify_original_extra_args(self, mock_getsize, mock_client_constructor):
        """Verify that the original extra_args dict is not modified (immutability)"""
        mock_client = MagicMock()
        mock_client_constructor.return_value = mock_client

        s3_provider = S3Provider()

        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            tmp_file.write(b"test data")
            tmp_file.flush()
            file_path = tmp_file.name

        try:
            mock_getsize.return_value = SIZE_50GB + 1

            original_extra_args = {"ACL": "bucket-owner-full-control"}
            original_extra_args_copy = original_extra_args.copy()

            s3_provider.upload_file(file_path, "test-bucket", "test-key", original_extra_args)

            # Verify original dict was not modified
            self.assertEqual(original_extra_args, original_extra_args_copy)
            self.assertNotIn("CacheControl", original_extra_args)

        finally:
            os.unlink(file_path)

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_size_50gb_constant(self, mock_client_constructor):
        """Verify SIZE_50GB constant is correct value"""
        expected_50gb = 50 * 1024 * 1024 * 1024  # 50GB in bytes
        self.assertEqual(SIZE_50GB, expected_50gb)
        self.assertEqual(SIZE_50GB, 53687091200)  # Explicit value check


if __name__ == "__main__":
    unittest.main()
