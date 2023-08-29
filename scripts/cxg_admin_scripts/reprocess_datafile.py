import json
import logging
import os
import sys
from time import time

import boto3
import click
from click import Context

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

logging.basicConfig()
logger: logging.Logger = logging.getLogger(__name__)


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
    """Cxg remaster v2"""
    pass


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
