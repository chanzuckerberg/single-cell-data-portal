import os
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import pytest

from backend.layers.thirdparty.s3_exceptions import S3DeleteException
from backend.layers.thirdparty.s3_provider import (
    CLOUDFRONT_MAX_CACHE_SIZE,
    S3Provider,
    _get_cloudfront_max_cache_size,
    _parse_size_string,
)


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
            # Mock file size to be > CloudFront max cache size
            mock_getsize.return_value = CLOUDFRONT_MAX_CACHE_SIZE + 1

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
            # Mock file size to be exactly at CloudFront max cache size threshold
            mock_getsize.return_value = CLOUDFRONT_MAX_CACHE_SIZE

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
            mock_getsize.return_value = CLOUDFRONT_MAX_CACHE_SIZE + 1

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
            mock_getsize.return_value = CLOUDFRONT_MAX_CACHE_SIZE + 1

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
            mock_getsize.return_value = CLOUDFRONT_MAX_CACHE_SIZE + 1

            original_extra_args = {"ACL": "bucket-owner-full-control"}
            original_extra_args_copy = original_extra_args.copy()

            s3_provider.upload_file(file_path, "test-bucket", "test-key", original_extra_args)

            # Verify original dict was not modified
            self.assertEqual(original_extra_args, original_extra_args_copy)
            self.assertNotIn("CacheControl", original_extra_args)

        finally:
            os.unlink(file_path)

    @patch("backend.layers.thirdparty.s3_provider.boto3.client")
    def test_cloudfront_max_cache_size_default(self, mock_client_constructor):
        """Verify CLOUDFRONT_MAX_CACHE_SIZE defaults to 50GB when env var not set"""
        expected_50gb = 50 * 1024 * 1024 * 1024  # 50GB in bytes
        # Default should be 50GB (53687091200 bytes)
        self.assertEqual(CLOUDFRONT_MAX_CACHE_SIZE, expected_50gb)
        self.assertEqual(CLOUDFRONT_MAX_CACHE_SIZE, 53687091200)  # Explicit value check


# Test size string parsing with parametrized tests
@pytest.mark.parametrize(
    "size_str,expected_bytes",
    [
        # GB tests
        ("50GB", 50 * 1024 * 1024 * 1024),
        ("1GB", 1024 * 1024 * 1024),
        ("100GB", 100 * 1024 * 1024 * 1024),
        ("50gb", 50 * 1024 * 1024 * 1024),  # Case insensitive
        ("50 Gb", 50 * 1024 * 1024 * 1024),  # With space
        # MB tests
        ("100MB", 100 * 1024 * 1024),
        ("1MB", 1024 * 1024),
        ("500mb", 500 * 1024 * 1024),  # Case insensitive
        # TB tests
        ("1TB", 1024**3),
        ("2TB", 2 * (1024**3)),
        ("1tb", 1024**3),  # Case insensitive
        # KB tests
        ("1024KB", 1024 * 1024),
        ("1KB", 1024),
        ("512kb", 512 * 1024),  # Case insensitive
        # Bytes tests
        ("1024B", 1024),
        ("1B", 1),
        ("100b", 100),  # Case insensitive
        # Decimal values
        ("1.5GB", int(1.5 * 1024 * 1024 * 1024)),
        ("0.5TB", int(0.5 * (1024**3))),
        # Whitespace handling
        (" 50GB ", 50 * 1024 * 1024 * 1024),
        ("50 GB", 50 * 1024 * 1024 * 1024),
        ("  100MB  ", 100 * 1024 * 1024),
    ],
)
def test_parse_size_string_valid(size_str, expected_bytes):
    """Test parsing valid size strings"""
    assert _parse_size_string(size_str) == expected_bytes


@pytest.mark.parametrize(
    "invalid_size_str",
    ["invalid", "50", "GB", "50XX", ""],
)
def test_parse_size_string_invalid(invalid_size_str):
    """Test that invalid size formats raise ValueError"""
    with pytest.raises(ValueError):
        _parse_size_string(invalid_size_str)


@pytest.mark.parametrize(
    "env_value,expected_bytes",
    [
        ("100GB", 100 * 1024 * 1024 * 1024),
        ("500MB", 500 * 1024 * 1024),
        ("1TB", 1024**3),
        ("50GB", 50 * 1024 * 1024 * 1024),
    ],
)
def test_parse_size_string_from_env_var(env_value, expected_bytes):
    """Test that size strings from environment variable are parsed correctly"""
    assert _parse_size_string(env_value) == expected_bytes


def test_cloudfront_max_cache_size_default():
    """Verify CLOUDFRONT_MAX_CACHE_SIZE defaults to 50GB when env var not set"""
    expected_50gb = 50 * 1024 * 1024 * 1024
    assert expected_50gb == CLOUDFRONT_MAX_CACHE_SIZE
    assert CLOUDFRONT_MAX_CACHE_SIZE == 53687091200  # Explicit value check


@patch("backend.layers.thirdparty.s3_provider.os.getenv")
def test_get_cloudfront_max_cache_size_with_env_var(mock_getenv):
    """Test that _get_cloudfront_max_cache_size uses environment variable when set"""
    mock_getenv.return_value = "100GB"
    result = _get_cloudfront_max_cache_size()
    assert result == 100 * 1024 * 1024 * 1024
    mock_getenv.assert_called_once_with("CLOUDFRONT_MAX_CACHE_SIZE")


@patch("backend.layers.thirdparty.s3_provider.os.getenv")
def test_get_cloudfront_max_cache_size_default_when_not_set(mock_getenv):
    """Test that _get_cloudfront_max_cache_size uses default when env var not set"""
    mock_getenv.return_value = None
    result = _get_cloudfront_max_cache_size()
    expected_50gb = 50 * 1024 * 1024 * 1024
    assert result == expected_50gb


@patch("backend.layers.thirdparty.s3_provider.os.getenv")
@patch("backend.layers.thirdparty.s3_provider.logger")
def test_get_cloudfront_max_cache_size_invalid_falls_back_to_default(mock_logger, mock_getenv):
    """Test that invalid env var format falls back to default 50GB"""
    mock_getenv.return_value = "invalid"
    result = _get_cloudfront_max_cache_size()
    expected_50gb = 50 * 1024 * 1024 * 1024
    assert result == expected_50gb
    mock_logger.warning.assert_called_once()


if __name__ == "__main__":
    unittest.main()
