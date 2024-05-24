import io
import logging
import os
import tempfile

import numpy as np
import pyvips
from PIL import Image

from backend.layers.thirdparty.s3_provider import S3Provider

logger = logging.getLogger(__name__)


class SpatialDataProcessor:
    """
    A class that processes spatial data and generates deep zoom assets.

    Args:
        s3_provider (S3Provider, optional): An instance of the S3Provider class. Defaults to None.

    Attributes:
        s3_provider (S3Provider): An instance of the S3Provider class.
        deployment_stage (str): The deployment stage (e.g., "staging").
        env (str): The environment based on the deployment stage.
        bucket_name (str): The name of the S3 bucket.
        asset_directory (str): The directory for the deep zoom assets.
    """

    def __init__(self, s3_provider: S3Provider = None):
        self.s3_provider = s3_provider if s3_provider else S3Provider()
        self.deployment_stage = os.getenv("DEPLOYMENT_STAGE", "staging")
        # relying on dev bucket for rdev and test
        self.env = "dev" if self.deployment_stage in ["rdev", "test"] else self.deployment_stage
        self.bucket_name = f"spatial-deepzoom-{self.env}"
        self.asset_directory = "spatial-deep-zoom"

    def _calculate_aspect_ratio_crop(self, image_size):
        """
        Calculate the aspect ratio crop coordinates for an image.

        Args:
            image_size (int, int): The size of the image (width, height).

        Returns:
            tuple: The crop coordinates (left, upper, right, lower).
        """
        width, height = image_size
        new_dimension = min(width, height)
        left = (width - new_dimension) / 2
        upper = (height - new_dimension) / 2
        right = (width + new_dimension) / 2
        lower = (height + new_dimension) / 2
        return tuple(map(int, (left, upper, right, lower)))

    def _fetch_image(self, content):
        """
        Prepare the image array for processing.

        Args:
            content (dict): The content dictionary containing the image array.

        Returns:
            np.ndarray: The prepared image array.
        """
        resolution = "fullres" if "fullres" in content["images"] else "hires"

        image_array = content["images"][resolution]
        image_array_uint8 = np.uint8(image_array)

        return image_array_uint8, resolution

    def _process_and_flip_image(self, image_array_uint8):
        """
        Process and flip the image.

        Args:
            image_array_uint8 (np.ndarray): The image array.

        Returns:
            np.ndarray: The processed and flipped image array.
        """
        Image.MAX_IMAGE_PIXELS = None  # Disable the image size limit
        try:
            with Image.fromarray(image_array_uint8) as img:
                # #####
                # Convert RGBA to RGB if needed for JPEG compatibility:
                # bake the alpha channel into the image by pasting it onto
                # a white background
                # #####
                if img.mode == "RGBA":
                    background = Image.new("RGB", img.size, (255, 255, 255))
                    background.paste(img, (0, 0), img)
                    img = background
                cropped_img = img.crop(self._calculate_aspect_ratio_crop(img.size))
                cropped_img.save(io.BytesIO(), format="JPEG", quality=100)

                # Flip the image vertically due to explorer client rendering images upside down
                flipped_img = cropped_img.transpose(Image.FLIP_TOP_BOTTOM)

                return np.array(flipped_img)
        except Exception:
            logger.exception("Error processing image")
            raise

    def _generate_deep_zoom_assets(self, image_array, assets_folder):
        """
        Generate deep zoom assets from the image array.

        Args:
            image_array (np.ndarray): The image array.
            assets_folder (str): The temporary directory to save the assets.
        """
        h, w, bands = image_array.shape
        linear = image_array.reshape(w * h * bands)
        image = pyvips.Image.new_from_memory(linear.data, w, h, bands, "uchar")
        image.dzsave(os.path.join(assets_folder, "spatial"), suffix=".jpeg")

    def _upload_assets(self, assets_folder):
        """
        Upload the deep zoom assets to the S3 bucket.

        Args:
            assets_folder (str): The folder containing the assets.
        """
        s3_uri = f"s3://{self.bucket_name}/{self.asset_directory}/{os.path.basename(assets_folder)}"
        self.s3_provider.upload_directory(assets_folder, s3_uri)

    def create_deep_zoom_assets(self, container_name, content):
        """
        Create deep zoom assets for a container.

        Args:
            container_name (str): The name of the container.
            content (dict): The content dictionary containing the image array.
        """
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                assets_folder = os.path.join(temp_dir, container_name.replace(".cxg", ""))
                os.makedirs(assets_folder)

                image_array, _ = self._fetch_image(content)
                processed_image = self._process_and_flip_image(image_array)
                self._generate_deep_zoom_assets(processed_image, assets_folder)
                self._upload_assets(assets_folder)
        except Exception as e:
            logger.exception(f"Failed to create and upload deep zoom assets: {e}")
            raise

    def filter_spatial_data(self, content, library_id):
        """
        Filter spatial data based on the library ID.

        Args:
            content (dict): The content dictionary.
            library_id (str): The library ID.

        Returns:
            dict: The filtered spatial data.
        """
        image_array_uint8, resolution = self._fetch_image(content)
        with Image.fromarray(image_array_uint8) as img:
            width, height = img.size
            width = height = min(width, height)
            crop_coords = self._calculate_aspect_ratio_crop(img.size)
        return {
            library_id: {
                "image_properties": {
                    "resolution": resolution,
                    "crop_coords": crop_coords,
                    "width": width,
                    "height": height,
                },
                "images": {"hires": content["images"]["hires"], "fullres": []},
                "scalefactors": {
                    "spot_diameter_fullres": content["scalefactors"]["spot_diameter_fullres"],
                    "tissue_hires_scalef": content["scalefactors"]["tissue_hires_scalef"],
                },
            }
        }
