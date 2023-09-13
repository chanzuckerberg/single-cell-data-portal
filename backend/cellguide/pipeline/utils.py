import gzip
import json
import logging
import os
import pathlib
from dataclasses import dataclass
from typing import Union

import pandas as pd

from backend.cellguide.pipeline.config import CellGuideConfig
from backend.cellguide.pipeline.constants import ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME
from backend.common.utils.dataclass import convert_dataclass_to_dict

logger = logging.getLogger(__name__)


def output_json(data, path):
    """
    Output JSON while handling the case where `data` is an instance of a dataclass

    Arguments
    ---------
    data - a dict or a dataclass
    path - str
    """

    data = convert_dataclass_to_dict_and_strip_nones(data)

    os.makedirs(os.path.dirname(path), exist_ok=True)

    if path.endswith(".gz"):
        with gzip.open(path, "wt") as f:
            json.dump(data, f)
    else:
        with open(path, "w") as f:
            json.dump(data, f)


def output_json_per_key(data, path):
    """
    Output each value in the input data to a separate JSON file.
    """

    data = convert_dataclass_to_dict_and_strip_nones(data)

    os.makedirs(os.path.dirname(f"{path}/"), exist_ok=True)

    for key in data:
        filename = key.replace(":", "_")
        output_json(data[key], f"{path}/{filename}.json")


def convert_dataclass_to_dict_and_strip_nones(data):
    """
    A wrapper for dataclass-to-dict conversion that deletes keys with None as values
    """
    return _remove_none_values(convert_dataclass_to_dict(data))


def _remove_none_values(data):
    """
    Recursively remove keys with None values.

    Arguments
    ---------
    data - dict
    """
    if isinstance(data, dict):
        return {k: _remove_none_values(v) for k, v in data.items() if v is not None}
    elif isinstance(data, list):
        return [_remove_none_values(i) for i in data]
    else:
        return data


@dataclass
class GeneMetadata:
    gene_id_to_name: dict[str, str]
    gene_id_to_symbol: dict[str, str]


def get_gene_id_to_name_and_symbol() -> GeneMetadata:
    """
    This function reads a file containing gene metadata and returns a GeneMetadata object.
    The GeneMetadata object contains two dictionaries: gene_id_to_name and gene_id_to_symbol.
    gene_id_to_name is a dictionary where the keys are Ensembl GeneIDs and the values are gene descriptions.
    gene_id_to_symbol is a dictionary where the keys are Ensembl GeneIDs and the values are gene symbols.

    Returns
    -------
    GeneMetadata
        An object containing two dictionaries: gene_id_to_name and gene_id_to_symbol.
    """

    file_path = pathlib.Path(__file__).absolute().parent.joinpath("fixtures", ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME)
    gene_metadata = pd.read_csv(file_path, sep="\t")
    gene_id_to_name = gene_metadata.set_index("Ensembl GeneIDs")["Description"].to_dict()
    gene_id_to_symbol = gene_metadata.set_index("Ensembl GeneIDs")["Symbols"].to_dict()
    return GeneMetadata(gene_id_to_name=gene_id_to_name, gene_id_to_symbol=gene_id_to_symbol)


def get_bucket_path() -> Union[str, None]:
    """
    This function generates a bucket path based on the deployment stage and remote development prefix.
    If the deployment stage is 'rdev' and a remote development prefix is provided, the bucket path is prefixed
    with 'env-rdev-cellguide' and the remote development prefix. Otherwise, the bucket path is the same as the
    bucket name.

    Returns
    -------
    str or None
        The generated bucket path if the deployment stage is set and valid, None otherwise.
    """

    deployment_stage = os.getenv("DEPLOYMENT_STAGE")
    remote_dev_prefix = os.getenv("REMOTE_DEV_PREFIX")

    if not deployment_stage:
        logger.warning("Not uploading the pipeline output to S3 as DEPLOYMENT_STAGE is not set.")
        return

    if deployment_stage == "rdev" and remote_dev_prefix is None:
        logger.warning(
            "Not uploading the pipeline output to S3 as REMOTE_DEV_PREFIX is not set when DEPLOYMENT_STAGE is rdev."
        )
        return
    elif deployment_stage == "rdev":
        bucket = CellGuideConfig().bucket
        bucket_path = f"s3://{bucket}/env-rdev-cellguide{remote_dev_prefix}/"
    elif deployment_stage in ["dev", "staging", "prod"]:
        bucket = CellGuideConfig().bucket
        bucket_path = f"s3://{bucket}/"
    else:
        logger.warning(
            f"Invalid DEPLOYMENT_STAGE value: {deployment_stage}. Please set DEPLOYMENT_STAGE to one of rdev, dev, staging, or prod"
        )
        return

    return bucket_path


def get_object_key(*, object: str) -> str:
    """
    This function generates an object key for S3 based on the deployment stage and remote development prefix.
    If the deployment stage is 'rdev' and a remote development prefix is provided, the object key is prefixed
    with 'env-rdev-cellguide' and the remote development prefix. Otherwise, the object key is the same as the
    input object.

    Parameters
    ----------
    object : str
        The input object for which the object key is to be generated.

    Returns
    -------
    str
        The generated object key.
    """

    deployment_stage = os.getenv("DEPLOYMENT_STAGE")
    remote_dev_prefix = os.getenv("REMOTE_DEV_PREFIX")

    object_key = object
    if deployment_stage == "rdev" and remote_dev_prefix is not None:
        CellGuideConfig().bucket
        object_key = f"env-rdev-cellguide{remote_dev_prefix}/{object}"

    return object_key
