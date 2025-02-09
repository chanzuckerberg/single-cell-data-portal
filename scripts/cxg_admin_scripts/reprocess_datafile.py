import logging
import os
import sys

import boto3

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

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
    """Cxg remaster v2"""
    pass
