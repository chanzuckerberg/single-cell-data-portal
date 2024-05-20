import os
import pickle
import tempfile
from os import mkdir
from uuid import uuid4

import numpy as np
import pytest
from PIL import Image

from backend.layers.processing.utils.cxg_generation_utils import convert_uns_to_cxg_group
from backend.layers.processing.utils.spatial import SpatialDataProcessor
from backend.layers.thirdparty.s3_provider import S3Provider
from tests.unit.backend.fixtures.environment_setup import fixture_file_path


@pytest.fixture
def testing_cxg_temp_directory():
    directory = fixture_file_path(str(uuid4()))
    mkdir(directory)
    yield directory
    # Cleanup can be added here if needed


@pytest.fixture
def s3_provider_mock(mocker):
    return mocker.MagicMock(spec=S3Provider)


@pytest.fixture
def valid_uns():
    return {
        "spatial": {
            "library_id_1": {
                "images": {"hires": np.random.rand(10, 10, 3), "fullres": np.random.rand(20, 20, 3)},
                "scalefactors": {"spot_diameter_fullres": 1.0, "tissue_hires_scalef": 0.5},
            }
        }
    }


@pytest.fixture
def valid_spatial_data(valid_uns):
    return valid_uns["spatial"]["library_id_1"]


@pytest.fixture
def invalid_uns():
    return {}


@pytest.fixture
def spatial_processor(s3_provider_mock):
    return SpatialDataProcessor(s3_provider=s3_provider_mock)


@pytest.fixture
def output_folder():
    return "test_output"


@pytest.fixture
def cxg_container():
    return "test_container"


@pytest.fixture
def group_metadata_name():
    return "uns"


@pytest.fixture
def mock_spatial_processor(mocker, valid_uns):
    processor = mocker.MagicMock(spec=SpatialDataProcessor)
    processor.filter_spatial_data.return_value = valid_uns["spatial"]
    return processor


@pytest.fixture
def library_id():
    return "abcd123"


@pytest.fixture
def ctx():
    return None


def test__valid_input_metadata_copy(spatial_processor, valid_spatial_data, library_id):
    """
    Test case for verifying the required metadata is present
    """
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


def test__invalid_input_metadata_copy(spatial_processor, invalid_uns, library_id):
    """
    Test case to verify that a KeyError is raised when passing invalid input metadata
    """
    with pytest.raises(KeyError):
        spatial_processor.filter_spatial_data(invalid_uns, library_id)


def test__fetch_image_fullres(spatial_processor, valid_spatial_data):
    """
    Test that _fetch_image returns the fullres image when present.
    """
    image_array = spatial_processor._fetch_image(valid_spatial_data)
    assert image_array.shape == (20, 20, 3), "Expected fullres image to be returned."


def test__fetch_image_hires(spatial_processor, valid_spatial_data):
    """
    Test that _fetch_image returns the hires image when fullres is not present.
    """
    valid_spatial_data_without_fullres = valid_spatial_data.copy()
    del valid_spatial_data_without_fullres["images"]["fullres"]

    image_array = spatial_processor._fetch_image(valid_spatial_data_without_fullres)
    assert image_array.shape == (10, 10, 3), "Expected hires image to be returned."


def test__fetch_image_key_error(spatial_processor, valid_spatial_data):
    """
    Test that _fetch_image raises KeyError when neither fullres nor hires images are present.
    """
    valid_spatial_data_without_images = valid_spatial_data.copy()
    del valid_spatial_data_without_images["images"]["fullres"]
    del valid_spatial_data_without_images["images"]["hires"]

    with pytest.raises(KeyError):
        spatial_processor._fetch_image(valid_spatial_data_without_images)


def test__process_and_flip_image(spatial_processor):
    """
    Test the image processing method
    """
    test_image_array = (np.random.rand(10, 10, 3) * 255).astype(np.uint8)
    expected_flipped_array = np.flipud(test_image_array)

    flipped_image_array = spatial_processor._process_and_flip_image(test_image_array)
    assert np.array_equal(flipped_image_array, expected_flipped_array)


def test__crop_to_aspect_ratio(spatial_processor):
    """
    Test the cropping method to ensure the image is cropped to a square
    """
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


def test__generate_deep_zoom_assets(spatial_processor, output_folder, mocker):
    """
    Test the method to generate deep zoom assets
    """
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


def test__upload_assets(spatial_processor, output_folder, mocker):
    """
    Test upload assets to S3
    """
    mock_upload = mocker.patch.object(spatial_processor.s3_provider, "upload_directory")

    spatial_processor._upload_assets(output_folder)
    expected_s3_uri = f"s3://{spatial_processor.bucket_name}/{spatial_processor.asset_directory}/{output_folder}"

    # verify that upload_directory was called correctly
    mock_upload.assert_called_once_with(output_folder, expected_s3_uri)


def test__upload_assets_failure(spatial_processor, output_folder, mocker):
    """
    Test upload assets to S3 when the upload fails
    """
    mock_upload = mocker.patch.object(spatial_processor.s3_provider, "upload_directory")
    mock_upload.side_effect = Exception("Failed to upload")

    with pytest.raises(Exception, match="Failed to upload"):
        spatial_processor._upload_assets(output_folder)

    expected_s3_uri = f"s3://{spatial_processor.bucket_name}/{spatial_processor.asset_directory}/{output_folder}"
    mock_upload.assert_called_once_with(output_folder, expected_s3_uri)


def test__create_deep_zoom_assets(spatial_processor, cxg_container, valid_spatial_data, mocker):
    mock_fetch_image = mocker.patch.object(spatial_processor, "_fetch_image")
    mock_process_and_flip_image = mocker.patch.object(spatial_processor, "_process_and_flip_image")
    mock_generate_deep_zoom_assets = mocker.patch.object(spatial_processor, "_generate_deep_zoom_assets")
    mock_upload_assets = mocker.patch.object(spatial_processor, "_upload_assets")

    # mock the TemporaryDirectory context manager
    mock_temp_dir = mocker.patch("tempfile.TemporaryDirectory")
    temp_dir_name = "/mock/temp/dir"
    mock_temp_dir.return_value.__enter__.return_value = temp_dir_name

    # mock return values for the internal methods
    mock_fetch_image.return_value = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)
    mock_process_and_flip_image.return_value = np.random.randint(0, 255, (100, 100, 3), dtype=np.uint8)

    # call the method under test
    spatial_processor.create_deep_zoom_assets(cxg_container, valid_spatial_data)

    assets_folder = os.path.join(temp_dir_name, cxg_container.replace(".cxg", ""))

    # assertions to ensure each step is called
    mock_fetch_image.assert_called_once_with(valid_spatial_data)
    mock_process_and_flip_image.assert_called_once_with(mock_fetch_image.return_value)
    mock_generate_deep_zoom_assets.assert_called_once_with(mock_process_and_flip_image.return_value, assets_folder)
    mock_upload_assets.assert_called_once_with(assets_folder)


def test__create_deep_zoom_assets_exception(spatial_processor, cxg_container, valid_spatial_data, mocker):
    mock_fetch_image = mocker.patch.object(spatial_processor, "_fetch_image")
    mock_process_and_flip_image = mocker.patch.object(spatial_processor, "_process_and_flip_image")
    mock_generate_deep_zoom_assets = mocker.patch.object(spatial_processor, "_generate_deep_zoom_assets")
    mock_upload_assets = mocker.patch.object(spatial_processor, "_upload_assets")

    # mock an exception in the _fetch_image method
    mock_fetch_image.side_effect = Exception("Test exception")

    # assert that the method raises an exception
    with pytest.raises(Exception, match="Test exception"):
        spatial_processor.create_deep_zoom_assets(cxg_container, valid_spatial_data)

    mock_fetch_image.assert_called_once_with(valid_spatial_data)
    mock_process_and_flip_image.assert_not_called()
    mock_generate_deep_zoom_assets.assert_not_called()
    mock_upload_assets.assert_not_called()


def test__convert_uns_to_cxg_group(cxg_container, valid_uns, group_metadata_name, ctx, mock_spatial_processor, mocker):
    mock_from_numpy = mocker.patch("backend.layers.processing.utils.cxg_generation_utils.tiledb.from_numpy")
    mock_tiledb_open = mocker.patch(
        "backend.layers.processing.utils.cxg_generation_utils.tiledb.open", mocker.mock_open()
    )
    mocker.patch(
        "backend.layers.processing.utils.cxg_generation_utils.SpatialDataProcessor",
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
