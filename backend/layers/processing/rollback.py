import os
from enum import Enum
from typing import List

from backend.layers.business.business import BusinessLogic

# from backend.layers.common.entities import (
# )
from backend.layers.persistence.persistence import DatabaseProvider
from backend.layers.thirdparty.s3_provider import S3Provider
from backend.layers.thirdparty.uri_provider import UriProvider


class RollbackType(Enum):
    PRIVATE_COLLECTIONS = "private_collections"
    PUBLIC_COLLECTIONS = "public_collections"
    COLLECTION_LIST = "collection_list"
    DATASET_LIST = "dataset_list"


class RollbackEntity:
    def __init__(
        self, business_logic: BusinessLogic, rollback_type: RollbackType, entity_id_list: List[str] = None
    ) -> None:
        self.business_logic = business_logic
        self.rollback_type = str(rollback_type)
        self.entity_id_list = entity_id_list

    def rollback(self):
        if self.rollback_type == RollbackType.PRIVATE_COLLECTIONS:
            self.rollback_private_collections()
        elif self.rollback_type == RollbackType.PUBLIC_COLLECTIONS:
            self.rollback_public_collections()
        elif self.rollback_type == RollbackType.COLLECTION_LIST:
            self.rollback_collection_list()
        elif self.rollback_type == RollbackType.DATASET_LIST:
            self.rollback_dataset_list()
        else:
            raise ValueError(f"Invalid rollback type: {self.rollback_type}")

    def rollback_dataset(self, dataset_version_id: str):
        pass

    def rollback_private_collection(self, collection_version_id: str):
        pass

    def rollback_public_collection(self, collection_id: str):
        pass

    def rollback_private_collections(self):
        pass

    def rollback_public_collections(self):
        pass

    def rollback_collection_list(self):
        pass

    def rollback_dataset_list(self):
        pass


if __name__ == "__main__":
    business_logic = BusinessLogic(
        DatabaseProvider(),
        None,
        None,
        None,
        S3Provider(),
        UriProvider(),
    )
    rollback_type = os.environ.get("ROLLBACK_TYPE", None)
    if rollback_type is None:
        raise ValueError("ROLLBACK_TYPE is required")
    if rollback_type not in RollbackType.__members__:
        raise ValueError(f"ROLLBACK_TYPE must be one of {', '.join(RollbackType.__members__)}")

    if rollback_type == RollbackType.COLLECTION_LIST:
        collections_to_rollback = os.environ.get("COLLECTION_LIST", None)
        RollbackEntity(business_logic, RollbackType(rollback_type), collections_to_rollback).rollback()
    elif rollback_type == RollbackType.DATASET_LIST:
        datasets_to_rollback = os.environ.get("DATASET_LIST", None)
        RollbackEntity(business_logic, RollbackType(rollback_type), datasets_to_rollback).rollback()
    else:
        RollbackEntity(business_logic, RollbackType(rollback_type)).rollback()
