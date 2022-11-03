from backend.layers.common.entities import CollectionVersionId, DatasetVersionId


class StepFunctionProviderInterface:
    def start_step_function(self, version_id: CollectionVersionId, dataset_version_id: DatasetVersionId, url: str) -> None:
        pass