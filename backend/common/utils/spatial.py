import io
import os
import shutil

import numpy as np
import pyvips
from PIL import Image

from backend.layers.thirdparty.s3_provider import S3Provider


class SpatialDataProcessor:
    def __init__(self, s3_provider: S3Provider):
        self.s3_provider = s3_provider if s3_provider else S3Provider()
        deployment_stage = os.getenv("DEPLOYMENT_STAGE", "staging")
        self.env = "dev" if deployment_stage in ["rdev", "test"] else deployment_stage
        self.bucket_name = f"spatial-deepzoom-{self.env}"
        self.asset_directory = "spatial-deep-zoom"
        self.base_dir = "./z_spatial/"

    def _calculate_aspect_ratio_crop(self, image_size):
        width, height = image_size
        new_dimension = min(width, height)
        left = (width - new_dimension) / 2
        upper = (height - new_dimension) / 2
        right = (width + new_dimension) / 2
        lower = (height + new_dimension) / 2
        # return (left, upper, left + new_dimension, upper + new_dimension)
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

    def _generate_deep_zoom_assets(self, image_array, container_name):
        h, w, bands = image_array.shape
        linear = image_array.reshape(w * h * bands)
        image = pyvips.Image.new_from_memory(linear.data, w, h, bands, "uchar")
        output_path = f"{self.base_dir}/{container_name}/"
        image.dzsave(output_path + "spatial", suffix=".jpeg")
        return output_path

    def _upload_assets(self, container_name, directory):
        s3_uri = f"s3://{self.bucket_name}/{self.asset_directory}/{container_name}/"
        self.s3_provider.upload_directory(directory, s3_uri)
        shutil.rmtree(directory)  # Cleanup the local directory after upload

    def create_deep_zoom_assets(self, container_name, content):
        image_array = self._prepare_image(content)
        processed_image = self._process_and_flip_image(image_array)
        output_dir = self._generate_deep_zoom_assets(processed_image, container_name)
        self._upload_assets(container_name, output_dir)

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
