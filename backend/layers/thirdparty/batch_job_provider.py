import json
import os
from time import time
from typing import Dict

import boto3

from backend.layers.common.entities import (
    DatasetArtifactMetadataUpdate,
    DatasetVersionId,
)


class BatchJobProviderInterface:
    def start_metadata_update_batch_job(
        self,
        current_dataset_version_id: DatasetVersionId,
        new_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ) -> None:
        pass


class BatchJobProvider(BatchJobProviderInterface):
    def __init__(self) -> None:
        self.client = boto3.client("batch")

    def start_metadata_update_batch_job(
        self,
        current_dataset_version_id: DatasetVersionId,
        new_dataset_version_id: DatasetVersionId,
        metadata_update: DatasetArtifactMetadataUpdate,
    ) -> Dict[str, str]:
        """
        Starts a Batch Job that updates metadata on all dataset artifacts with indicated mapped changes
        """
        deployment_stage = os.environ.get("DEPLOYMENT_STAGE")
        stack_name = "stagestack" if deployment_stage == "staging" else f"{deployment_stage}stack"
        if os.environ.get("REMOTE_DEV_PREFIX"):
            stack_name = os.environ.get("REMOTE_DEV_PREFIX").replace("/", "")
        return self.client.submit_job(
            jobName=f"metadata_update_{new_dataset_version_id}_{int(time())}",
            jobQueue=f"dp-{deployment_stage}",
            jobDefinition=f"dp-{deployment_stage}-{stack_name}-dataset-metadata-update",
            containerOverrides={
                "environment": [
                    {
                        "name": "CURRENT_DATASET_VERSION_ID",
                        "value": current_dataset_version_id.id,
                    },
                    {
                        "name": "NEW_DATASET_VERSION_ID",
                        "value": new_dataset_version_id.id,
                    },
                    {
                        "name": "METADATA_UPDATE_JSON",
                        "value": json.dumps(metadata_update.as_dict_without_none_values()),
                    },
                ]
            },
        )
