import copy
import os
import pickle
import tempfile
from os import mkdir
from uuid import uuid4

import numpy as np
import pytest
from PIL import Image

from backend.common.utils.cxg_generation_utils import convert_uns_to_cxg_group
from backend.layers.processing.spatial import SpatialDataProcessor
from backend.layers.thirdparty.s3_provider import S3Provider
from tests.unit.backend.fixtures.environment_setup import fixture_file_path


@pytest.fixture
def setup_spatial_processor(mocker):
    testing_cxg_temp_directory = fixture_file_path(str(uuid4()))

    s3_provider_mock = mocker.MagicMock(spec=S3Provider)
    library_id = "abcd123"
    valid_uns = {
        "spatial": {
            "library_id_1": {
                "images": {"hires": np.random.rand(10, 10, 3), "fullres": np.random.rand(20, 20, 3)},
                "scalefactors": {"spot_diameter_fullres": 1.0, "tissue_hires_scalef": 0.5},
            }
        }
    }
    valid_spatial_data = valid_uns["spatial"]["library_id_1"]
    invalid_uns = {}
    spatial_processor = SpatialDataProcessor(s3_provider=s3_provider_mock)
    mkdir(testing_cxg_temp_directory)

    output_folder = "test_output"
    cxg_container = "test_container"
    group_metadata_name = "uns"
    ctx = None

    mock_spatial_processor = mocker.MagicMock(spec=SpatialDataProcessor)
    mock_spatial_processor.filter_spatial_data.return_value = valid_uns["spatial"]

    return {
        "testing_cxg_temp_directory": testing_cxg_temp_directory,
        "s3_provider_mock": s3_provider_mock,
        "library_id": library_id,
        "valid_uns": valid_uns,
        "valid_spatial_data": valid_spatial_data,
        "invalid_uns": invalid_uns,
        "spatial_processor": spatial_processor,
        "output_folder": output_folder,
        "cxg_container": cxg_container,
        "group_metadata_name": group_metadata_name,
        "ctx": ctx,
        "mock_spatial_processor": mock_spatial_processor,
    }


def test__valid_input_metadata_copy(setup_spatial_processor):
    """
    Test case for verifying the required metadata is present
    """
    spatial_processor = setup_spatial_processor["spatial_processor"]
    valid_spatial_data = setup_spatial_processor["valid_spatial_data"]
    library_id = setup_spatial_processor["library_id"]

    # Expected output based on the input data
    expected_output = {
        library_id: {
            "images": {"hires": valid_spatial_data["images"]["hires"], "fullres": []},
            "scalefactors": {
                "spot_diameter_fullres": valid_spatial_data["scalefactors"]["spot_diameter_fullres"],
                "tissue_hires_scalef": valid_spatial_data["scalefactors"]["tissue_hires_scalef"],
            },
        }
    }

    result = spatial_processor.filter_spatial_data(valid_spatial_data, library_id)

    assert result == expected_output, f"Expected {expected_output}, but got {result}"


def test__invalid_input_metadata_copy(setup_spatial_processor):
    """
    Test case to verify that a KeyError is raised when passing input metadata
    """
    spatial_processor = setup_spatial_processor["spatial_processor"]
    invalid_uns = setup_spatial_processor["invalid_uns"]
    library_id = setup_spatial_processor["library_id"]
    with pytest.raises(KeyError):
        spatial_processor.filter_spatial_data(invalid_uns, library_id)


def test__full_high_res_image_present(setup_spatial_processor):
    """
    Test that the default image is fullres if present otherwise fall back to hires
    """
    spatial_processor = setup_spatial_processor["spatial_processor"]
    valid_spatial_data = setup_spatial_processor["valid_spatial_data"]

    image_array_uint8 = spatial_processor._prepare_image(valid_spatial_data)
    assert image_array_uint8.shape == (20, 20, 3)

    valid_hires_uns = copy.deepcopy(valid_spatial_data)
    valid_hires_uns["images"].pop("fullres", None)
    image_array_uint8 = spatial_processor._prepare_image(valid_hires_uns)
    assert image_array_uint8.shape == (10, 10, 3)

    non_valid_uns = copy.deepcopy(valid_hires_uns)
    non_valid_uns["images"].pop("hires", None)
    with pytest.raises(KeyError):
        image_array_uint8 = spatial_processor._prepare_image(non_valid_uns)


def test__process_and_flip_image(setup_spatial_processor):
    """
    Test the image processing method
    """
    spatial_processor = setup_spatial_processor["spatial_processor"]
    test_image_array = (np.random.rand(10, 10, 3) * 255).astype(np.uint8)
    expected_flipped_array = np.flipud(test_image_array)

    flipped_image_array = spatial_processor._process_and_flip_image(test_image_array)
    assert np.array_equal(flipped_image_array, expected_flipped_array)


def test__crop_to_aspect_ratio(setup_spatial_processor):
    """
    Test the cropping method to ensure the image is cropped to a square
    """
    spatial_processor = setup_spatial_processor["spatial_processor"]
    width, height = 23, 11
    test_image_array = (np.random.rand(height, width, 3) * 255).astype(np.uint8)

    test_image = Image.fromarray(test_image_array)

    crop_box = spatial_processor._calculate_aspect_ratio_crop(test_image.size)
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


def test__generate_deep_zoom_assets(setup_spatial_processor, mocker):
    """
    Test the method to generate deep zoom assets
    """
    spatial_processor = setup_spatial_processor["spatial_processor"]
    output_folder = setup_spatial_processor["output_folder"]
    test_image_array = (np.random.rand(100, 100, 3) * 255).astype(np.uint8)

    mock_new_from_memory = mocker.patch("pyvips.Image.new_from_memory")
    mock_image = mocker.MagicMock()
    mock_new_from_memory.return_value = mock_image

    with tempfile.TemporaryDirectory() as temp_dir:
        assets_folder = os.path.join(temp_dir, output_folder)
        os.makedirs(assets_folder)
        spatial_processor._generate_deep_zoom_assets(test_image_array, assets_folder)

        # verify new_from_memory was called correctly
        h, w, bands = test_image_array.shape
        linear = test_image_array.reshape(w * h * bands)
        mock_new_from_memory.assert_called_once_with(linear.data, w, h, bands, "uchar")

        # verify dzsave was called correctly on the mock_image object
        expected_output_path = os.path.join(assets_folder, "spatial")
        mock_image.dzsave.assert_called_once_with(expected_output_path, suffix=".webp")


def test__upload_assets(setup_spatial_processor, mocker):
    """
    Test upload assets to S3
    """
    spatial_processor = setup_spatial_processor["spatial_processor"]
    output_folder = setup_spatial_processor["output_folder"]

    mock_upload = mocker.patch.object(spatial_processor.s3_provider, "upload_directory")

    spatial_processor._upload_assets(output_folder)
    expected_s3_uri = f"s3://{spatial_processor.bucket_name}/{spatial_processor.asset_directory}/{output_folder}"

    # verify that upload_directory was called correctly
    mock_upload.assert_called_once_with(output_folder, expected_s3_uri)


def test__upload_assets_failure(setup_spatial_processor, mocker):
    """
    Test upload assets to S3 when the upload fails
    """
    spatial_processor = setup_spatial_processor["spatial_processor"]
    output_folder = setup_spatial_processor["output_folder"]

    mock_upload = mocker.patch.object(spatial_processor.s3_provider, "upload_directory")
    mock_upload.side_effect = Exception("Failed to upload")

    with pytest.raises(Exception, match="Failed to upload"):
        spatial_processor._upload_assets(output_folder)

    expected_s3_uri = f"s3://{spatial_processor.bucket_name}/{spatial_processor.asset_directory}/{output_folder}"
    mock_upload.assert_called_once_with(output_folder, expected_s3_uri)


def test__create_deep_zoom_assets(setup_spatial_processor, mocker):
    spatial_processor = setup_spatial_processor["spatial_processor"]
    cxg_container = setup_spatial_processor["cxg_container"]
    valid_spatial_data = setup_spatial_processor["valid_spatial_data"]

    mock_prepare_image = mocker.patch.object(spatial_processor, "_prepare_image")
    mock_process_and_flip_image = mocker.patch.object(spatial_processor, "_process_and_flip_image")
    mock_generate_deep_zoom_assets = mocker.patch.object(spatial_processor, "_generate_deep_zoom_assets")
    mock_upload_assets = mocker.patch.object(spatial_processor, "_upload_assets")

    # mock the TemporaryDirectory context manager
    mock_temp_dir = mocker.patch("tempfile.TemporaryDirectory")
    temp_dir_name = "/mock/temp/dir"
    mock_temp_dir.return_value.__enter__.return_value = temp_dir_name

    # mock return values for the internal methods
    mock_prepare_image.return_value = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)
    mock_process_and_flip_image.return_value = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)

    # call the method under test
    spatial_processor.create_deep_zoom_assets(cxg_container, valid_spatial_data)

    assets_folder = os.path.join(temp_dir_name, cxg_container.replace(".cxg", ""))

    # assertions to ensure each step is called
    mock_prepare_image.assert_called_once_with(valid_spatial_data)
    mock_process_and_flip_image.assert_called_once_with(mock_prepare_image.return_value)
    mock_generate_deep_zoom_assets.assert_called_once_with(mock_process_and_flip_image.return_value, assets_folder)
    mock_upload_assets.assert_called_once_with(assets_folder)


def test__create_deep_zoom_assets_exception(setup_spatial_processor, mocker):
    spatial_processor = setup_spatial_processor["spatial_processor"]
    cxg_container = setup_spatial_processor["cxg_container"]
    valid_spatial_data = setup_spatial_processor["valid_spatial_data"]

    mock_prepare_image = mocker.patch.object(spatial_processor, "_prepare_image")
    mock_process_and_flip_image = mocker.patch.object(spatial_processor, "_process_and_flip_image")
    mock_generate_deep_zoom_assets = mocker.patch.object(spatial_processor, "_generate_deep_zoom_assets")
    mock_upload_assets = mocker.patch.object(spatial_processor, "_upload_assets")

    # mock an exception in the _prepare_image method
    mock_prepare_image.side_effect = Exception("Test exception")

    # assert that the method raises an exception
    with pytest.raises(Exception, match="Test exception"):
        spatial_processor.create_deep_zoom_assets(cxg_container, valid_spatial_data)

    mock_prepare_image.assert_called_once_with(valid_spatial_data)
    mock_process_and_flip_image.assert_not_called()
    mock_generate_deep_zoom_assets.assert_not_called()
    mock_upload_assets.assert_not_called()


def test__convert_uns_to_cxg_group(setup_spatial_processor, mocker):
    cxg_container = setup_spatial_processor["cxg_container"]
    valid_uns = setup_spatial_processor["valid_uns"]
    group_metadata_name = setup_spatial_processor["group_metadata_name"]
    ctx = setup_spatial_processor["ctx"]
    mock_spatial_processor = setup_spatial_processor["mock_spatial_processor"]

    mock_from_numpy = mocker.patch("backend.common.utils.cxg_generation_utils.tiledb.from_numpy")
    mock_tiledb_open = mocker.patch("backend.common.utils.cxg_generation_utils.tiledb.open", mocker.mock_open())
    mocker.patch(
        "backend.common.utils.cxg_generation_utils.SpatialDataProcessor",
        return_value=mock_spatial_processor,
    )

    mock_metadata_array = mock_tiledb_open.return_value.__enter__.return_value
    mock_metadata_array.meta = {}

    convert_uns_to_cxg_group(cxg_container, valid_uns, group_metadata_name, ctx)

    # check if from_numpy is called correctly
    mock_from_numpy.assert_called_once_with(f"{cxg_container}/{group_metadata_name}", np.zeros((1,)))

    # check if spatial metadata is processed and written correctly
    assert "spatial" in mock_metadata_array.meta
    spatial_data = pickle.loads(mock_metadata_array.meta["spatial"])
    assert "library_id_1" in spatial_data
    assert np.array_equal(
        spatial_data["library_id_1"]["images"]["hires"],
        valid_uns["spatial"]["library_id_1"]["images"]["hires"],
    )
    assert (
        spatial_data["library_id_1"]["scalefactors"]["spot_diameter_fullres"]
        == valid_uns["spatial"]["library_id_1"]["scalefactors"]["spot_diameter_fullres"]
    )
    assert (
        spatial_data["library_id_1"]["scalefactors"]["tissue_hires_scalef"]
        == valid_uns["spatial"]["library_id_1"]["scalefactors"]["tissue_hires_scalef"]
    )

    # check if create_deep_zoom_assets is called correctly
    mock_spatial_processor.create_deep_zoom_assets.assert_called_once_with(
        cxg_container, valid_uns["spatial"]["library_id_1"]
    )
