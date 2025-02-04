import json
import os

import jsonschema

local_path = os.path.dirname(os.path.realpath(__file__))
with open(local_path + "/schema.json") as fp:
    schema = json.load(fp)
validator = jsonschema.Draft7Validator(schema)


def to_manifest(anndata: str, atac_seq_fragment: str = None) -> dict:
    manifest = {"anndata": anndata}
    if atac_seq_fragment:
        manifest["atac_seq_fragment"] = atac_seq_fragment
    return manifest
