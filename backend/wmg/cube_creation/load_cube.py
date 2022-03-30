import logging
import os
import sys
import time
from typing import Dict
import subprocess
from pathlib import Path

import tiledb

import pandas as pd

from backend.corpora.common.corpora_orm import DatasetArtifactFileType
from backend.corpora.common.entities import Dataset, DatasetAsset
from backend.corpora.common.utils.db_session import db_session_manager
from backend.wmg.cube_creation.corpus_schema import create_tdb
from backend.wmg.cube_creation.loader import load
from backend.wmg.cube_creation.wmg_cube import create_cube

from pronto import Ontology
import pygraphviz as pgv
import json

import boto3

# TODO - does writing and reading directly from s3 slow down compute? test

Log_Format = "%(levelname)s %(asctime)s - %(message)s"

logging.basicConfig(stream=sys.stdout, filemode="w", format=Log_Format, level=logging.WARNING)

logger = logging.getLogger(__name__)

# Stack name for rdev
stack_name = os.environ.get("REMOTE_DEV_PREFIX")
wmg_bucket_name = os.environ.get("WMG_BUCKET")
artifact_bucket_name = os.environ.get("ARTIFACT_BUCKET")


def get_wmg_bucket_path():
    if stack_name:
        return f"s3://{wmg_bucket_name}{stack_name}"
    else:
        return f"s3://{wmg_bucket_name}"


def get_s3_uris():
    with db_session_manager() as session:
        dataset_ids = Dataset.list_ids_for_cube(session)
        s3_uris = DatasetAsset.s3_uris_for_datasets(session, dataset_ids, DatasetArtifactFileType.H5AD)
    return s3_uris


def copy_datasets_to_instance(dataset_directory):
    s3_uris = get_s3_uris()
    for dataset in s3_uris.keys():
        sync_command = ["aws", "s3", "cp", s3_uris[dataset], f"./{dataset_directory}/{dataset}/local.h5ad"]
    subprocess.run(sync_command)  # TODO parallelize this step


def load_datasets_into_corpus(path_to_datasets, group_name):
    # try:
    load(path_to_datasets, group_name, True)
    # except Exception as e:
    #     logger.error(f"Issue loading datasets into corpus: {e}")


def get_cells_by_tissue_type(tdb_group: str) -> Dict:
    with tiledb.open(f"{tdb_group}/obs", "r") as obs:
        cell_tissue_types = (
            obs.query(attrs=[], dims=["tissue_ontology_term_id", "cell_type_ontology_term_id"])
            .df[:]
            .drop_duplicates()
            .sort_values(by="tissue_ontology_term_id")
        )
    unique_tissue_ontology_term_id = cell_tissue_types.tissue_ontology_term_id.unique()
    cell_type_by_tissue = {}
    for x in unique_tissue_ontology_term_id:
        cell_type_by_tissue[x] = cell_tissue_types.loc[
            cell_tissue_types["tissue_ontology_term_id"] == x, "cell_type_ontology_term_id"
        ]

    return cell_type_by_tissue


def generate_cell_ordering(cell_type_by_tissue):
    onto = Ontology.from_obo_library("cl-basic.obo")

    def compute_ordering(cells, root):
        ancestors = [list(onto[t].superclasses()) for t in cells if t in onto]
        ancestors = [i for s in ancestors for i in s]
        ancestors = set(ancestors)

        G = pgv.AGraph()
        for a in ancestors:
            for s in a.subclasses(with_self=False, distance=1):
                if s in ancestors:
                    G.add_edge(a.id, s.id)

        G.layout(prog="dot")

        positions = {}
        for n in G.iternodes():
            pos = n.attr["pos"].split(",")
            positions[n] = (float(pos[0]), float(pos[1]))

        ancestor_ids = [a.id for a in ancestors]

        def recurse(node):
            if node in cells:
                yield (node)
            children = [
                (c, positions[c.id]) for c in onto[node].subclasses(with_self=False, distance=1) if c.id in ancestor_ids
            ]
            sorted_children = sorted(children, key=lambda x: x[1][0])
            for child in sorted_children:
                yield from recurse(child[0].id)

        ordered_list = list(dict.fromkeys(recurse(root)))
        return ordered_list

    mapping = {}
    for tissue, cell_df in cell_type_by_tissue.items():
        cells = list(cell_df)
        ordered_cells = compute_ordering(cells, "CL:0000003")  # TODO: is this the right root?
        mapping[tissue] = ordered_cells

    data = []
    for tissue, cells in mapping.items():
        for cell in cells:
            data.append((tissue, cell))

    df = pd.DataFrame(data, columns=["tissue", "cell"])
    df.to_json("ordered-cells.json")


def update_s3_resources(group_name):
    timestamp = int(time.time())
    upload_cube_to_s3(group_name, timestamp)
    upload_cell_ordering_to_s3(timestamp)
    update_latest_snapshot_identifier(timestamp)
    remove_oldest_datasets(timestamp)


def remove_oldest_datasets(timestamp):
    s3 = boto3.resource("s3")
    wmg_bucket = s3.Bucket(wmg_bucket_name)
    objects = wmg_bucket.objects.filter(Prefix=stack_name.strip("/")) if stack_name else wmg_bucket.objects.all()

    def enumerate_timestamps_and_objects(objects):
        for obj in objects:
            try:
                tokens = obj.key.split("/")
                timestamp_token = tokens[1] if stack_name else tokens[0]
                is_timestamp = timestamp_token[:10].isdigit()
                if is_timestamp:
                    yield (timestamp_token, obj)
            except Exception:
                pass

    candidate_to_delete = [obj for obj in enumerate_timestamps_and_objects(objects)]

    timestamps = sorted(list(set([x[0] for x in candidate_to_delete])))

    if len(timestamps) > 2:
        timestamps_to_delete = list(timestamps)[:-2]
    else:
        timestamps_to_delete = []

    for timestamp, object in candidate_to_delete:
        if timestamp in timestamps_to_delete:
            object.delete()


def upload_cube_to_s3(group_name, timestamp):
    sync_command = ["aws", "s3", "sync", f"{group_name}/cube", f"{get_wmg_bucket_path()}/{timestamp}/cube"]
    subprocess.run(sync_command)


def upload_cell_ordering_to_s3(timestamp):
    sync_command = ["aws", "s3", "cp", f"ordered-cells.json", f"{get_wmg_bucket_path()}/{timestamp}/ordered-cells.json"]
    subprocess.run(sync_command)


def update_latest_snapshot_identifier(timestamp):
    s3 = boto3.resource("s3")
    file_key = "latest_snapshot_identifier"
    file_path = f"{stack_name.strip('/')}/{file_key}" if stack_name else f"{file_key}"
    object = s3.Object(wmg_bucket_name, file_path)
    object.put(Body=str(timestamp))


def load_data_and_create_cube(path_to_datasets, group_name):
    if not tiledb.VFS().is_dir(group_name):
        create_tdb(group_name)
    copy_datasets_to_instance(path_to_datasets)
    load_datasets_into_corpus(path_to_datasets, group_name)
    try:
        create_cube(group_name)
    except Exception as e:
        logger.error(f"Issue creating the cube: {e}")
    cell_type_by_tissue = get_cells_by_tissue_type(group_name)
    generate_cell_ordering(cell_type_by_tissue)
    update_s3_resources(group_name)
    print("Cube creation script - completed")
    return 0


if __name__ == "__main__":
    rv = load_data_and_create_cube("datasets", "full-cube")
    sys.exit(rv)
