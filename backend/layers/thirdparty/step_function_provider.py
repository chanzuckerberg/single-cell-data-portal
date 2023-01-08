import json
from time import time

import boto3

from backend.common.corpora_config import CorporaConfig
from backend.layers.common.entities import CollectionVersionId, DatasetVersionId


class StepFunctionProviderInterface:
    def start_step_function(
        self, version_id: CollectionVersionId, dataset_version_id: DatasetVersionId, url: str
    ) -> None:
        pass


class StepFunctionProvider(StepFunctionProviderInterface):
    def __init__(self) -> None:
        self.client = boto3.client("stepfunctions")

    def start_step_function(
        self, version_id: CollectionVersionId, dataset_version_id: DatasetVersionId, url: str
    ) -> None:
        input_parameters = {
            "collection_id": version_id.id,
            "url": url,
            "dataset_id": dataset_version_id.id,
        }
        sfn_name = f"{dataset_version_id}_{int(time())}"
        response = self.client.start_execution(
            stateMachineArn=CorporaConfig().upload_sfn_arn,
            name=sfn_name,
            input=json.dumps(input_parameters),
        )
        return response
