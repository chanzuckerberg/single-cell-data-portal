import boto3
import click
from click import Context

import logging
import json
import os
import sys
from time import time, sleep

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.corpora_orm import (
    DbDataset,
    DatasetArtifactFileType,
)
from backend.corpora.common.entities.dataset import Dataset

from urllib.parse import urlparse

logging.basicConfig()
logger = logging.getLogger(__name__)


def get_aws_account_id() -> str:
    sts = boto3.client("sts")
    return sts.get_caller_identity()["Account"]


def get_happy_stack_name(deployment) -> str:
    """
    Returns the name of the Happy stack for the specified deployment
    Note: This will only work with deployment={dev,stage,prod} and will not work with rdev!
    :param deployment: dev, stage or prod
    :return:
    """
    return f"{deployment}-{deployment}stack"


def cxg_remaster(ctx):
    """Cxg remaster"""

    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to remaster all the cxgs?",
            abort=True,
        )

        client = boto3.client("stepfunctions")

        for record in session.query(DbDataset):
            if not record.tombstone:
                dataset = Dataset.get(session, dataset_id=record.id)
                artifacts = [a.s3_uri for a in dataset.artifacts if a.filetype == DatasetArtifactFileType.CXG]
                if len(artifacts) > 0:
                    cxg = artifacts[0]

                    p = urlparse(cxg)
                    bucket = p.hostname
                    dataset_id = p.path.strip("/").strip(".cxg")

                    if dataset_id == "2e5273bd-aa36-4478-8f6f-62fa0abcea43":
                        continue

                    if bucket != "hosted-cellxgene-dev":
                        continue

                    print(bucket, dataset_id)

                    input = {"dataset_id": dataset_id}

                    aws_account_id = get_aws_account_id()
                    deployment = ctx.obj["deployment"]
                    happy_stack_name = get_happy_stack_name(deployment)

                    response = client.start_execution(
                        stateMachineArn=f"arn:aws:states:us-west-2:{aws_account_id}:stateMachine:dp-"
                        f"{happy_stack_name}-cxg-remaster-sfn",
                        name=f"{dataset_id}-{int(time())}",
                        input=json.dumps(input),
                    )

                    print(response["executionArn"])
                    sleep(1)

def cxg_remaster_v2(ctx):
    """Cxg remaster v2"""

    with db_session_manager() as session:
        click.confirm(
            "Are you sure you want to remaster all the cxgs?",
            abort=True,
        )

        client = boto3.client("stepfunctions")

        for record in session.query(DbDataset):
            if not record.tombstone:
                dataset = Dataset.get(session, dataset_id=record.id)
                artifacts = [a.s3_uri for a in dataset.artifacts if a.filetype == DatasetArtifactFileType.CXG]
                if len(artifacts) > 0:
                    cxg = artifacts[0]

                    p = urlparse(cxg)
                    bucket = p.hostname
                    dataset_id = p.path.strip("/").strip(".cxg")

                    if dataset_id == "2e5273bd-aa36-4478-8f6f-62fa0abcea43":
                        continue

                    if bucket != "hosted-cellxgene-dev":
                        continue

                    print(bucket, dataset_id)

                    input = {"dataset_id": dataset_id}

                    aws_account_id = get_aws_account_id()
                    deployment = ctx.obj["deployment"]
                    happy_stack_name = get_happy_stack_name(deployment)

                    response = client.start_execution(
                        stateMachineArn=f"arn:aws:states:us-west-2:{aws_account_id}:stateMachine:dp-"
                        f"{happy_stack_name}-cxg-remaster-v2-sfn",
                        name=f"{dataset_id}-{int(time())}",
                        input=json.dumps(input),
                    )

                    print(response["executionArn"])
                    sleep(1)


def reprocess_seurat(ctx: Context, dataset_id: str) -> None:
    """
    Reconverts the specified dataset to Seurat format in place.
    :param ctx: command context
    :param dataset_id: ID of dataset to reconvert to Seurat format
    """

    deployment = ctx.obj["deployment"]

    click.confirm(
        f"Are you sure you want to run this script? "
        f"It will reconvert and replace the dataset {dataset_id} to Seurat in the {deployment} environment.",
        abort=True,
    )

    aws_account_id = get_aws_account_id()
    deployment = ctx.obj["deployment"]
    happy_stack_name = get_happy_stack_name(deployment)

    payload = {"dataset_id": dataset_id}

    client = boto3.client("stepfunctions")
    response = client.start_execution(
        stateMachineArn=f"arn:aws:states:us-west-2:{aws_account_id}:stateMachine:dp-{happy_stack_name}-seurat-sfn",
        name=f"{dataset_id}-{int(time())}",
        input=json.dumps(payload),
    )

    click.echo(
        f"Step function executing: "
        f"https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/executions/details/"
        f"{response['executionArn']}"
    )
