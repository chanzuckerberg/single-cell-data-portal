# TODO(mbarrien): This should not be calling TestDatabase unconditionally upon import.

from ..backend.corpora.fixtures.database import TestDatabase


TestDatabase()
