#! /usr/bin/env python

import os
import sys
from time import time

import boto3
import click

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))  # noqa
sys.path.insert(0, pkg_root)  # noqa


@click.group()
def cli():
    pass


@cli.command()
@click.argument("deployment_stage", type=click.Choice(["dev", "staging", "prod", "rdev"]))
@click.argument(
    "rollback_type",
    type=click.Choice(
        [
            "private_dataset_list",
            "private_collection_list",
            "private_collections",
            "public_collection_list",
            "public_collections",
        ]
    ),
)
@click.option(
    "--entity-id-list",
    default="",
    help="comma-delimited ID list of entities to rollback, if using '*_list' rollback type",
)
@click.option("--rdev-stack-name", default="", help="if deployment stage is rdev, specific stack name")
def trigger_rollback(deployment_stage, rollback_type, entity_id_list, rdev_stack_name):
    """
    Used to rollback datasets or collections to a previous version.

    ./scripts/cxg_admin_scripts/rollback.py trigger-rollback <deployment_stage> <rollback_type> --entity-id-list entity_id1,entity_id2
    """
    rollback(rollback_type, deployment_stage, entity_id_list, rdev_stack_name)


def rollback(rollback_type: str, deployment_stage: str, entity_id_list: str = "", rdev_stack_name: str = ""):
    stack_name = f"{deployment_stage}stack"
    if deployment_stage == "staging":
        stack_name = "stagestack"
    elif deployment_stage == "rdev":
        stack_name = rdev_stack_name

    client = boto3.client("batch")
    response = client.submit_job(
        jobName=f"rollback_{int(time())}",
        jobQueue=f"schema_migration-{deployment_stage}",
        jobDefinition=f"dp-{deployment_stage}-{stack_name}-rollback",
        containerOverrides={
            "environment": [
                {
                    "name": "ROLLBACK_TYPE",
                    "value": rollback_type,
                },
                {
                    "name": "ENTITY_LIST",
                    "value": entity_id_list,
                },
            ]
        },
    )

    click.echo(
        f"Batch Job executing: "
        f"https://us-west-2.console.aws.amazon.com/states/home?region=us-west-2#/executions/details/"
        f"{response['jobArn']}"
    )


if __name__ == "__main__":
    cli()
