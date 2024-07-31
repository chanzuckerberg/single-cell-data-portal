import json
import os
from typing import List, TypedDict

import boto3

from backend.layers.common.entities import DatasetVersionId

AsyncLambdaInvocationResponse = TypedDict("AsyncLambdaInvocationResponse", {"StatusCode": int})


class LambdaProviderInterface:
    def invoke_dataset_version_cleanup_handler(
        self, dataset_version_ids: List[DatasetVersionId]
    ) -> AsyncLambdaInvocationResponse:  # type: ignore
        pass


class LambdaProvider(LambdaProviderInterface):
    def __init__(self) -> None:
        self.client = boto3.client("lambda")

    def invoke_dataset_version_cleanup_handler(
        self, dataset_version_ids: List[DatasetVersionId]
    ) -> AsyncLambdaInvocationResponse:
        """
        Starts a lambda that will ingest a list of DatasetVersion ids and delete those
        DatasetVersions and associated artifacts in the db and S3
        """
        deployment_stage = os.environ.get("DEPLOYMENT_STAGE")
        stack_name = "stagestack" if deployment_stage == "staging" else f"{deployment_stage}stack"
        if os.environ.get("REMOTE_DEV_PREFIX") is not None:
            stack_name = os.environ.get("REMOTE_DEV_PREFIX").replace("/", "")
        function_name = f"dp-{deployment_stage}-{stack_name}-dataset-version-cleanup"

        lambda_payload = {
            "dataset_version_ids": [str(version_id) for version_id in dataset_version_ids],
        }
        response: AsyncLambdaInvocationResponse = self.client.invoke(
            FunctionName=function_name, InvocationType="Event", Payload=json.dumps(lambda_payload)
        )
        return response
