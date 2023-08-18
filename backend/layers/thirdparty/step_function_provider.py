import json
from time import time

import boto3

from backend.common.corpora_config import CorporaConfig
from backend.layers.common.entities import CollectionVersionId, DatasetVersionId


def sfn_name_generator(dataset_version_id: DatasetVersionId, prefix=None) -> str:
    if prefix:
        return f"{prefix}_{dataset_version_id}_{int(time())}"
    else:
        return f"{dataset_version_id}_{int(time())}"


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
        """
        Starts a step function that will ingest the dataset `dataset_version_id` using the artifact
        located at `url`
        """
        input_parameters = {
            "collection_id": version_id.id,
            "url": url,
            "dataset_id": dataset_version_id.id,
        }
        sfn_name = sfn_name_generator(dataset_version_id)
        response = self.client.start_execution(
            stateMachineArn=CorporaConfig().upload_sfn_arn,
            name=sfn_name,
            input=json.dumps(input_parameters),
        )
        return response
