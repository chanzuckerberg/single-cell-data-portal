from http import HTTPStatus

from server.common.errors import DatasetAccessError
from server.common.utils.data_locator import DataLocator


class DataLoader(object):
    def __init__(self, location, app_config=None):
        """location can be a string or DataLocator"""
        self.app_config = app_config
        self.dataset_config = self.__resolve_dataset_config()
        region_name = None if app_config is None else app_config.server__data_locator__s3_region_name
        self.location = DataLocator(location, region_name=region_name)
        if not self.location.exists():
            raise DatasetAccessError("Dataset does not exist.", HTTPStatus.NOT_FOUND)

        from server.dataset.cxg_dataset import CxgDataset

        self.matrix_type = CxgDataset

    def __resolve_dataset_config(self):
        dataset_config = self.app_config.default_dataset
        if dataset_config is None:
            raise DatasetAccessError("Missing dataset config", HTTPStatus.NOT_FOUND)
        return dataset_config

    def pre_load_validation(self):
        self.matrix_type.pre_load_validation(self.location)

    def file_size(self):
        return self.matrix_type.file_size(self.location)

    def open(self):
        # create and return a DataAdaptor object
        return self.matrix_type.open(self.location, self.app_config)

    def validate_and_open(self):
        # create and return a DataAdaptor object
        self.pre_load_validation()
        return self.open()
