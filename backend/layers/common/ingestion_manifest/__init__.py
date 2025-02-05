import json
import os

import jsonschema


def get_schema() -> dict:
    local_path = os.path.dirname(os.path.realpath(__file__))
    with open(local_path + "/schema.json") as fp:
        return json.load(fp)


def get_validator() -> jsonschema.Draft202012Validator:
    schema = get_schema()
    return jsonschema.Draft202012Validator(schema)


def to_manifest(anndata: str, atac_seq_fragment: str = None) -> dict:
    manifest = {"anndata": anndata}
    if atac_seq_fragment:
        manifest["atac_seq_fragment"] = atac_seq_fragment
    return manifest
