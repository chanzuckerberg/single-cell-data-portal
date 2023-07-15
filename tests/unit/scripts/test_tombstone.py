from unittest.mock import MagicMock

from click import Context

from scripts.cxg_admin_scripts.tombstones import tombstone_collection
from tests.unit.backend.layers.common.base_test import BaseTest


class TestTombstone(BaseTest):
    """
    Test scripting for tombstoning Collections and Datasets
    """

    def test__tombstone_collection(self):
        c_v = self.generate_published_collection()
        self.assertFalse(c_v.canonical_collection.tombstoned)
        [self.assertFalse(d.canonical_dataset.tombstoned) for d in c_v.datasets]

        context = MagicMock(spec=Context)
        context.obj = {"business_logic": self.business_logic}
        tombstone_collection(context, c_v.canonical_collection.id.id)

        c_v = self.business_logic.get_collection_version(c_v.version_id, get_tombstoned=True)
        self.assertTrue(c_v.canonical_collection.tombstoned)
