import boto3

from backend.cellguide.common.constants import CELL_GUIDE_VALID_EXPLORER_CXGS_FILENAME
from backend.cellguide.pipeline.explorer_cxgs.constants import (
    CXG_BUCKET_NAME,
    CXG_ORGANISMS,
    CXG_PATH_PREFIX,
    CXG_TISSUE_PATH_PREFIX,
)
from backend.cellguide.pipeline.explorer_cxgs.types import ValidExplorerCxgs
from backend.cellguide.pipeline.utils import output_json


def run(*, output_directory: str):
    valid_explorer_cxgs = get_valid_cxgs()
    output_json(valid_explorer_cxgs, f"{output_directory}/{CELL_GUIDE_VALID_EXPLORER_CXGS_FILENAME}")


def get_folders_from_s3(bucket, prefix):
    s3_client = boto3.client("s3")

    continuation_token = None

    all_objects = []
    while True:
        # Check if continuation token is present
        if continuation_token:
            response = s3_client.list_objects_v2(
                Bucket=bucket, Prefix=prefix, Delimiter="/", ContinuationToken=continuation_token
            )
        else:
            response = s3_client.list_objects_v2(Bucket=bucket, Prefix=prefix, Delimiter="/")
        # Extend the list of all objects with the current page of results
        if "CommonPrefixes" in response:
            all_objects.extend(response["CommonPrefixes"])

        # Check if there are more pages of results
        if response.get("IsTruncated"):
            continuation_token = response.get("NextContinuationToken")
        else:
            break

    return all_objects


def get_valid_cxgs() -> ValidExplorerCxgs:

    organism_cxgs = {}

    for organism in CXG_ORGANISMS:
        results = get_folders_from_s3(CXG_BUCKET_NAME, CXG_PATH_PREFIX + organism + "/")
        cxgs = [i["Prefix"].split(".cxg")[0].split("/")[-1].replace("_", ":") for i in results if "CL_" in i["Prefix"]]
        organism_cxgs[organism] = cxgs

    organism_tissue_cxgs = {}

    for organism in CXG_ORGANISMS:
        organism_tissue_cxgs[organism] = {}
        results = get_folders_from_s3(CXG_BUCKET_NAME, CXG_TISSUE_PATH_PREFIX + organism + "/")
        cxgs = [i["Prefix"].split(".cxg")[0].split("/")[-1] for i in results]
        cxgs = [[j.replace("_", ":") for j in i.split("__")] for i in cxgs if len(i.split("__")) == 2]
        for i in cxgs:
            t, c = i
            L = organism_tissue_cxgs[organism].get(t, set())
            L.add(c)
            organism_tissue_cxgs[organism][t] = L

        organism_tissue_cxgs[organism] = {k: list(v) for k, v in organism_tissue_cxgs[organism].items()}

    return ValidExplorerCxgs(organism_celltype_cxgs=organism_cxgs, organism_tissue_celltype_cxgs=organism_tissue_cxgs)
