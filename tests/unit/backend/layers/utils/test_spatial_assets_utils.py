import copy
import pickle
import unittest
from os import mkdir
from unittest.mock import MagicMock, mock_open, patch
from uuid import uuid4

import numpy as np
import pytest
from PIL import Image

from backend.common.utils.cxg_generation_utils import convert_uns_to_cxg_group
from backend.common.utils.spatial import SpatialDataProcessor
from backend.layers.thirdparty.s3_provider import S3Provider
from tests.unit.backend.fixtures.environment_setup import fixture_file_path


class TestSpatialDataProcessor(unittest.TestCase):
    def setUp(self):
        self.testing_cxg_temp_directory = fixture_file_path(str(uuid4()))

        self.s3_provider_mock = MagicMock(spec=S3Provider)
        self.library_id = "abcd123"
        self.valid_uns = {
            "spatial": {
                "library_id_1": {
                    "images": {"hires": np.random.rand(10, 10, 3), "fullres": np.random.rand(20, 20, 3)},
                    "scalefactors": {"spot_diameter_fullres": 1.0, "tissue_hires_scalef": 0.5},
                }
            }
        }
        self.valid_spatial_data = self.valid_uns["spatial"]["library_id_1"]
        self.invalid_uns = {}
        self.spatial_processor = SpatialDataProcessor(s3_provider=self.s3_provider_mock)
        mkdir(self.testing_cxg_temp_directory)

        self.output_folder = "test_output"
        self.cxg_container = "test_container"
        self.group_metadata_name = "uns"
        self.ctx = None

        self.mock_spatial_processor = MagicMock(spec=SpatialDataProcessor)
        self.mock_spatial_processor.filter_spatial_data.return_value = {
            "library_id_1": {
                "images": {"hires": self.valid_uns["spatial"]["library_id_1"]["images"]["hires"], "fullres": []},
                "scalefactors": {
                    "spot_diameter_fullres": self.valid_uns["spatial"]["library_id_1"]["scalefactors"][
                        "spot_diameter_fullres"
                    ],
                    "tissue_hires_scalef": self.valid_uns["spatial"]["library_id_1"]["scalefactors"][
                        "tissue_hires_scalef"
                    ],
                },
            }
        }

    def test__valid_input_metadata_copy(self):
        """
        Test case for verifying the required metadata is present
        """
        self.spatial_processor.filter_spatial_data(self.valid_spatial_data, self.library_id)

    def test__invalid_input_metadata_copy(self):
        """
        Test case to verify that a KeyError is raised when passing input metadata
        """
        with pytest.raises(KeyError):
            self.spatial_processor.filter_spatial_data(self.invalid_uns, self.library_id)

    def test__full_high_res_image_present(self):
        """
        Test that the default image is fullres if present otherwsie fall back to hires
        """
        image_array_uint8 = self.spatial_processor._prepare_image(self.valid_spatial_data)
        self.assertEqual(image_array_uint8.shape, (20, 20, 3))

        valid_hires_uns = copy.deepcopy(self.valid_spatial_data)
        valid_hires_uns["images"].pop("fullres", None)
        image_array_uint8 = self.spatial_processor._prepare_image(valid_hires_uns)
        self.assertEqual(image_array_uint8.shape, (10, 10, 3))

        non_valid_uns = copy.deepcopy(valid_hires_uns)
        non_valid_uns["images"].pop("hires", None)
        with pytest.raises(KeyError):
            image_array_uint8 = self.spatial_processor._prepare_image(non_valid_uns)

    def test__process_and_flip_image(self):
        """
        Test the image processing method
        """
        test_image_array = (np.random.rand(10, 10, 3) * 255).astype(np.uint8)
        expected_flipped_array = np.flipud(test_image_array)

        flipped_image_array = self.spatial_processor._process_and_flip_image(test_image_array)
        assert np.array_equal(flipped_image_array, expected_flipped_array)

    def test__crop_to_aspect_ratio(self):
        """
        Test the cropping method to ensure the image is cropped to a square
        """
        width, height = 20, 10
        test_image_array = (np.random.rand(height, width, 3) * 255).astype(np.uint8)

        test_image = Image.fromarray(test_image_array)

        crop_box = self.spatial_processor._calculate_aspect_ratio_crop(test_image.size)
        cropped_image = test_image.crop(crop_box)

        assert cropped_image.width == cropped_image.height, "Cropped image is not square"
        assert cropped_image.width == min(test_image.size), "Cropped dimension is incorrect"

        # verify that the crop is centered
        expected_width = expected_height = min(test_image.size)
        expected_left = (test_image.width - expected_width) / 2
        expected_upper = (test_image.height - expected_height) / 2
        assert crop_box == (
            expected_left,
            expected_upper,
            expected_left + expected_width,
            expected_upper + expected_height,
        ), "Crop is not centered correctly"

    def test__generate_deep_zoom_assets(self):
        """
        Test the method to generate deep zoom assets
        """
        test_image_array = (np.random.rand(100, 100, 3) * 255).astype(np.uint8)

        with patch("pyvips.Image.new_from_memory") as mock_new_from_memory:
            mock_image = MagicMock()
            mock_new_from_memory.return_value = mock_image

            self.spatial_processor._generate_deep_zoom_assets(test_image_array, self.output_folder)

            # verify new_from_memory was called correctly
            h, w, bands = test_image_array.shape
            linear = test_image_array.reshape(w * h * bands)
            mock_new_from_memory.assert_called_once_with(linear.data, w, h, bands, "uchar")

            # verify dzsave was called correctly on the mock_image object
            expected_output_path = self.output_folder + "spatial"
            mock_image.dzsave.assert_called_once_with(expected_output_path, suffix=".jpeg")

    def test__upload_assets(self):
        """
        Test upload assets to S3
        """

        with (
            patch.object(self.spatial_processor.s3_provider, "upload_directory") as mock_upload,
            patch("shutil.rmtree") as mock_rmtree,
        ):

            self.spatial_processor._upload_assets(self.output_folder)
            expected_s3_uri = f"s3://{self.spatial_processor.bucket_name}/{self.spatial_processor.asset_directory}/{self.output_folder}"

            # verify that upload_directory was called correctly
            mock_upload.assert_called_once_with(self.output_folder, expected_s3_uri)

            # verify that shutil.rmtree was called to remove the local directory
            mock_rmtree.assert_called_once_with(self.output_folder)

    def test__upload_assets_failure(self):
        """
        Test upload assets to S3 when the upload fails
        """

        with (
            patch.object(self.spatial_processor.s3_provider, "upload_directory") as mock_upload,
            patch("shutil.rmtree") as mock_rmtree,
        ):
            mock_upload.side_effect = Exception("Failed to upload")

            with pytest.raises(Exception) as exc_info:
                self.spatial_processor._upload_assets(self.output_folder)

            assert str(exc_info.value) == "Failed to upload"

            expected_s3_uri = f"s3://{self.spatial_processor.bucket_name}/{self.spatial_processor.asset_directory}/{self.output_folder}"
            mock_upload.assert_called_once_with(self.output_folder, expected_s3_uri)
            mock_rmtree.assert_not_called()  # assuming the directory is not deleted if upload fails

    def test__create_deep_zoom_assets(self):
        with (
            patch.object(self.spatial_processor, "_prepare_image") as mock_prepare_image,
            patch.object(self.spatial_processor, "_process_and_flip_image") as mock_process_and_flip_image,
            patch.object(self.spatial_processor, "_generate_deep_zoom_assets") as mock_generate_deep_zoom_assets,
            patch.object(self.spatial_processor, "_upload_assets") as mock_upload_assets,
        ):
            # Mock return values for the internal methods
            mock_prepare_image.return_value = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)
            mock_process_and_flip_image.return_value = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)

            # Call the method under test
            self.spatial_processor.create_deep_zoom_assets(self.cxg_container, self.valid_spatial_data)

            # Assertions to ensure each step is called
            mock_prepare_image.assert_called_once_with(self.valid_spatial_data)
            mock_process_and_flip_image.assert_called_once_with(mock_prepare_image.return_value)
            mock_generate_deep_zoom_assets.assert_called_once_with(
                mock_process_and_flip_image.return_value, "test_container/"
            )
            mock_upload_assets.assert_called_once_with("test_container/")

    def test__create_deep_zoom_assets_exception(self):
        with (
            patch.object(self.spatial_processor, "_prepare_image") as mock_prepare_image,
            patch.object(self.spatial_processor, "_process_and_flip_image") as mock_process_and_flip_image,
            patch.object(self.spatial_processor, "_generate_deep_zoom_assets") as mock_generate_deep_zoom_assets,
            patch.object(self.spatial_processor, "_upload_assets") as mock_upload_assets,
        ):
            # Mock an exception in the _prepare_image method
            mock_prepare_image.side_effect = Exception("Test exception")

            # Assert that the method raises an exception
            with self.assertRaises(Exception) as context:
                self.spatial_processor.create_deep_zoom_assets(self.cxg_container, self.valid_spatial_data)

            self.assertIn("An error occurred while creating and uploading deep zoom assets", str(context.exception))
            mock_prepare_image.assert_called_once_with(self.valid_spatial_data)
            mock_process_and_flip_image.assert_not_called()
            mock_generate_deep_zoom_assets.assert_not_called()
            mock_upload_assets.assert_not_called()

    def test__convert_uns_to_cxg_group(self):
        with (
            patch("backend.common.utils.cxg_generation_utils.tiledb.from_numpy") as mock_from_numpy,
            patch("backend.common.utils.cxg_generation_utils.tiledb.open", mock_open()) as mock_tiledb_open,
            patch(
                "backend.common.utils.cxg_generation_utils.SpatialDataProcessor",
                return_value=self.mock_spatial_processor,
            ),
        ):

            mock_metadata_array = mock_tiledb_open.return_value.__enter__.return_value
            mock_metadata_array.meta = {}

            convert_uns_to_cxg_group(self.cxg_container, self.valid_uns, self.group_metadata_name, self.ctx)

            # Check if from_numpy is called correctly
            mock_from_numpy.assert_called_once_with(f"{self.cxg_container}/{self.group_metadata_name}", np.zeros((1,)))

            # Check if metadata is processed and written correctly
            self.assertIn("spatial", mock_metadata_array.meta)
            spatial_data = pickle.loads(mock_metadata_array.meta["spatial"])
            self.assertIn("library_id_1", spatial_data)
            self.assertEqual(
                spatial_data["library_id_1"]["images"]["hires"].all(),
                self.valid_uns["spatial"]["library_id_1"]["images"]["hires"].all(),
            )
            self.assertEqual(
                spatial_data["library_id_1"]["scalefactors"]["spot_diameter_fullres"],
                self.valid_uns["spatial"]["library_id_1"]["scalefactors"]["spot_diameter_fullres"],
            )
            self.assertEqual(
                spatial_data["library_id_1"]["scalefactors"]["tissue_hires_scalef"],
                self.valid_uns["spatial"]["library_id_1"]["scalefactors"]["tissue_hires_scalef"],
            )

            # Check if create_deep_zoom_assets is called correctly
            self.mock_spatial_processor.create_deep_zoom_assets.assert_called_once_with(
                self.cxg_container, self.valid_uns["spatial"]["library_id_1"]
            )
