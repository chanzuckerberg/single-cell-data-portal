import json
import os
from time import time
from typing import Dict

import boto3

from backend.layers.common.entities import (
    CollectionVersionId,
    DatasetVersionId,
)


class BatchJobProviderInterface:
    def start_metadata_update_batch_job(
        self,
        collection_version_id: CollectionVersionId,
        dataset_version_id: DatasetVersionId,
        metadata_update_dict: Dict[str, str],
    ) -> None:
        pass


class BatchJobProvider(BatchJobProviderInterface):
    def __init__(self) -> None:
        self.client = boto3.client("batch")

    def start_metadata_update_batch_job(
        self,
        collection_version_id: CollectionVersionId,
        dataset_version_id: DatasetVersionId,
        metadata_update_dict: Dict[str, str],
    ) -> Dict[str, str]:
        """
        Starts a Batch Job that updates metadata on all dataset artifacts with indicated mapped changes
        """
        deployment_stage = os.environ.get("DEPLOYMENT_STAGE")
        stack_name = "stagestack" if deployment_stage == "staging" else f"{deployment_stage}stack"
        if os.environ.get("REMOTE_DEV_PREFIX"):
            stack_name = os.environ.get("REMOTE_DEV_PREFIX").replace("/", "")
        return self.client.submit_job(
            jobName=f"metadata_update_{dataset_version_id}_{int(time())}",
            jobQueue=f"dp-{deployment_stage}",
            jobDefinition=f"dp-{deployment_stage}-{stack_name}-dataset-metadata-update",
            containerOverrides={
                "environment": [
                    {
                        "name": "COLLECTION_VERSION_ID",
                        "value": collection_version_id.id,
                    },
                    {
                        "name": "DATASET_VERSION_ID",
                        "value": dataset_version_id.id,
                    },
                    {
                        "name": "METADATA_UPDATE_JSON",
                        "value": json.dumps(metadata_update_dict),
                    },
                ]
            },
        )
