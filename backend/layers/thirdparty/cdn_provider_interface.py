from typing import Protocol


class CDNProviderInterface(Protocol):
    def create_invalidation_for_index_paths(self):
        pass
