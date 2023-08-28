from unittest.mock import MagicMock

import pytest
from click import Context

from scripts.cxg_admin_scripts.tombstones import resurrect_collection, tombstone_collection
from tests.unit.backend.layers.common.base_test import BaseTest


class TestTombstone(BaseTest):
    """
    Test scripting for tombstoning Collections and Datasets
    """

    @pytest.fixture
    def context(self):
        context = MagicMock(spec=Context)
        context.obj = {"business_logic": self.business_logic}
        return context

    def test__tombstone_collection(self, context):
        c_v = self.generate_published_collection()

        tombstone_collection(context, c_v.collection_id.id)

        c_v = self.business_logic.get_collection_version(c_v.version_id, get_tombstoned=True)
        self.assertTrue(c_v.canonical_collection.tombstoned)

    def test__resurrect_collection(self, context):
        c_v = self.generate_published_collection()

        tombstone_collection(context, c_v.canonical_collection.id.id)

        c_v = self.business_logic.get_collection_version(c_v.version_id, get_tombstoned=True)

        resurrect_collection(context, c_v.collection_id.id)
        c_v = self.business_logic.get_collection_version(c_v.version_id)
        self.assertFalse(c_v.canonical_collection.tombstoned)
