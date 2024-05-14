import io
import os
import shutil

import numpy as np
import pyvips
from PIL import Image

from backend.layers.thirdparty.s3_provider import S3Provider


class SpatialDataProcessor:
    def __init__(self, s3_provider: S3Provider = None):
        self.s3_provider = s3_provider if s3_provider else S3Provider()
        self.deployment_stage = os.getenv("DEPLOYMENT_STAGE", "staging")
        # relying on dev bucket for rdev and test
        self.env = "dev" if self.deployment_stage in ["rdev", "test"] else self.deployment_stage
        self.bucket_name = f"spatial-deepzoom-{self.env}"
        self.asset_directory = "spatial-deep-zoom"

    def _calculate_aspect_ratio_crop(self, image_size):
        width, height = image_size
        new_dimension = min(width, height)
        left = (width - new_dimension) / 2
        upper = (height - new_dimension) / 2
        right = (width + new_dimension) / 2
        lower = (height + new_dimension) / 2
        return (left, upper, right, lower)

    def _prepare_image(self, content):
        image_key = "fullres" if "fullres" in content["images"] else "hires"
        image_array = content["images"][image_key]
        image_array_uint8 = np.uint8(image_array * 255 if image_key == "hires" else image_array)
        return image_array_uint8

    def _process_and_flip_image(self, image_array_uint8):
        Image.MAX_IMAGE_PIXELS = None  # Disable the image size limit
        with Image.fromarray(image_array_uint8) as img:
            cropped_img = img.crop(self._calculate_aspect_ratio_crop(img.size))  # Crop the image
            cropped_img.save(io.BytesIO(), format="JPEG", quality=90)  # Save or manipulate as needed
            flipped_img = cropped_img.transpose(Image.FLIP_TOP_BOTTOM)  # Flip the image vertically
        return np.array(flipped_img)

    def _generate_deep_zoom_assets(self, image_array, folder_name):
        h, w, bands = image_array.shape
        linear = image_array.reshape(w * h * bands)
        image = pyvips.Image.new_from_memory(linear.data, w, h, bands, "uchar")
        image.dzsave(folder_name + "spatial", suffix=".jpeg")

    def _upload_assets(self, assets_folder):
        s3_uri = f"s3://{self.bucket_name}/{self.asset_directory}/{assets_folder}"
        self.s3_provider.upload_directory(assets_folder, s3_uri)
        shutil.rmtree(assets_folder)

    def create_deep_zoom_assets(self, container_name, content):
        assets_folder = container_name.replace(".cxg", "") + "/"
        image_array = self._prepare_image(content)
        processed_image = self._process_and_flip_image(image_array)
        self._generate_deep_zoom_assets(processed_image, assets_folder)
        self._upload_assets(assets_folder)

    def filter_spatial_data(self, content, library_id):
        return {
            library_id: {
                "images": {"hires": content["images"]["hires"], "fullres": []},
                "scalefactors": {
                    "spot_diameter_fullres": content["scalefactors"]["spot_diameter_fullres"],
                    "tissue_hires_scalef": content["scalefactors"]["tissue_hires_scalef"],
                },
            }
        }
