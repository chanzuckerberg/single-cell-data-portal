"""Unit tests to verify guidescan2 is installed in the processing docker image."""

import shutil
import subprocess
import unittest


class TestGuidescanInstallation(unittest.TestCase):
    """Tests to verify guidescan2 is installed and accessible."""

    def test_guidescan_in_path(self):
        """Test that guidescan binary is available in PATH."""
        guidescan_path = shutil.which("guidescan")
        self.assertIsNotNone(
            guidescan_path, "guidescan binary not found in PATH. Ensure it is installed in the Docker image."
        )

    def test_guidescan_help_command(self):
        """Test that guidescan --help command executes successfully."""
        result = subprocess.run(
            ["guidescan", "--help"],
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(
            result.returncode,
            0,
            f"guidescan --help failed with return code {result.returncode}. "
            f"stderr: {result.stderr}, stdout: {result.stdout}",
        )

        # Verify the output contains expected content
        self.assertIn("Guidescan", result.stdout, "guidescan help output should contain 'Guidescan'")
        self.assertIn("Usage:", result.stdout, "guidescan help output should contain 'Usage:'")

    def test_guidescan_version_command(self):
        """Test that guidescan --version command executes successfully."""
        result = subprocess.run(
            ["guidescan", "--version"],
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(
            result.returncode,
            0,
            f"guidescan --version failed with return code {result.returncode}. "
            f"stderr: {result.stderr}, stdout: {result.stdout}",
        )

        # Verify version output contains version information
        self.assertGreater(
            len(result.stdout.strip()),
            0,
            "guidescan --version should output version information",
        )
