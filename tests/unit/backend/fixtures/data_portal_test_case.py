from backend.layers.common.entities import CollectionVersion
from tests.unit.backend.layers.common.base_test import BaseTest


class DataPortalTestCase(BaseTest):
    def generate_collection(self, **params) -> CollectionVersion:  # type: ignore
        visibility = params.pop("visibility", "PUBLIC")
        if visibility == "PUBLIC":
            return self.generate_published_collection(**params)  # type: ignore
        else:
            return self.generate_unpublished_collection(**params)
