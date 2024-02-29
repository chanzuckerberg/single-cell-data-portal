import boto3

from backend.cellguide.pipeline.constants import CELL_GUIDE_VALID_EXPLORER_CXGS_FILENAME
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


def get_valid_cxgs() -> ValidExplorerCxgs:
    s3_client = boto3.client("s3")

    organism_cxgs = {}

    for organism in CXG_ORGANISMS:
        results = s3_client.list_objects_v2(
            Bucket=CXG_BUCKET_NAME, Prefix=CXG_PATH_PREFIX + organism + "/", Delimiter="/"
        )
        cxgs = [i["Prefix"].split(".cxg")[0].split("/")[-1].replace("_", ":") for i in results["CommonPrefixes"]]
        organism_cxgs[organism] = cxgs

    organism_tissue_cxgs = {}

    for organism in CXG_ORGANISMS:
        organism_tissue_cxgs[organism] = {}
        results = s3_client.list_objects_v2(
            Bucket=CXG_BUCKET_NAME, Prefix=CXG_TISSUE_PATH_PREFIX + organism + "/", Delimiter="/"
        )
        cxgs = [i["Prefix"].split(".cxg")[0].split("/")[-1] for i in results["CommonPrefixes"]]
        cxgs = [[j.replace("_", ":") for j in i.split("__")] for i in cxgs]
        for i in cxgs:
            t, c = i
            L = organism_tissue_cxgs[organism].get(t, set())
            L.add(c)
            organism_tissue_cxgs[organism][t] = L

        organism_tissue_cxgs[organism] = {k: list(v) for k, v in organism_tissue_cxgs[organism].items()}

    return ValidExplorerCxgs(organism_celltype_cxgs=organism_cxgs, organism_tissue_celltype_cxgs=organism_tissue_cxgs)
