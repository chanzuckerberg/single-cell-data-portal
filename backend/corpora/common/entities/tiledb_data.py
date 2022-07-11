import os
import shutil
import numpy as np
import tiledb
import uuid

"""
SCHEMA WITH TILEDB VERSIONING:

s3://single-cell-corpus/

    metadata/ [GROUP]
        collections [ARRAY]
            - uuid
            - owner
            - published
            - revision_of
            - datasets
            - other metadata
        datasets [ARRAY]
            - uuid
            - name
            - other metadata

    datasets/<UUID>/
        artifacts/
        <UUID>.h5ad
        <UUID>.rds
        <UUID>.cxg
    wmg/<SNAPSHOT>/
    soma/<SNAPSHOT>/
"""

def new_id():
    return uuid.uuid4().hex

class TileDBData():
    # for testing purposes
    @staticmethod
    def init_db(location):
        # create group
        if os.path.exists(location):
            shutil.rmtree(location)

        tiledb.group_create(location)

        # create collections array
        # TODO: figure out ideal domain and tile
        # TODO: should we use more than one dimension?
        dim1 = tiledb.Dim(name="uuid", domain=(None, None), tile=2, dtype="S0")
        dom = tiledb.Domain(dim1)
        # TODO: add rest of required attributes
        a1 = tiledb.Attr(name="owner", dtype="U1")
        a2 = tiledb.Attr(name="published", dtype=np.int32) # 0=private, 1=published, -1=deleted
        a3 = tiledb.Attr(name="name", dtype="U1")
        a4 = tiledb.Attr(name="description", dtype="U1")
        a5 = tiledb.Attr(name="contact_name", dtype="U1")
        a6 = tiledb.Attr(name="contact_email", dtype="U1")
        a7 = tiledb.Attr(name="links", dtype="U1")
        a8 = tiledb.Attr(name="datasets", dtype="U1")
        a9 = tiledb.Attr(name="revision_of", dtype="U1")

        schema = tiledb.ArraySchema(domain=dom, sparse=True, attrs=[a1, a2, a3, a4, a5, a6, a7, a8, a9])
        array = location + "/collections"
        tiledb.Array.create(array, schema)

        # create datasets array
        # TODO: add more attributes as needed
        a1 = tiledb.Attr(name="name", dtype="U1")
        a2 = tiledb.Attr(name="artifact_id", dtype="U1")
        schema = tiledb.ArraySchema(domain=dom, sparse=True, attrs=[a1, a2])
        array = location + "/datasets"
        tiledb.Array.create(array, schema)

    # for testing purposes
    @staticmethod
    def destroy_db(location):
        shutil.rmtree(location)

    def __init__(self, location):
        self.location = location

    def create_collection(
        self,
        name: str = "",
        description: str = "",
        owner: str = "",
        contact_name: str = "",
        contact_email: str = "",
        links: list = None,
        datasets: list = None
    ):
        id = new_id()
        with tiledb.open(self.location + "/collections", "w") as A:
            A[id] = {
                "published": False,
                "name": name,
                "description": description,
                "contact_name": contact_name,
                "contact_email": contact_email,
                "links": ",".join(links), # TODO: figure out actual list support in TileDB
                "owner": owner,
                "datasets": ",".join(datasets),
                "revision_of": ""
            }
        return id

    def get_collection(self, id):
        with tiledb.open(self.location + "/collections", 'r') as A:
            return A[id]

    def get_attribute(self, id, attr):
        coll = self.get_collection(id)
        data = coll[attr][0]
        if attr == "links" or attr == "datasets":
            data = data.split(',') if len(data) > 0 else []
        return data

    def edit_collection(self, id, key, val):
        new_data = None
        with tiledb.open(self.location + "/collections", "r") as A:
            data = A[id]
            new_data = {
                "published": data['published'][0],
                "name": data['name'][0],
                "description": data['description'][0],
                "contact_name": data['contact_name'][0],
                "contact_email": data['contact_email'][0],
                "links": data['links'][0], # TODO: figure out actual list support in TileDB
                "owner": data['owner'][0],
                "datasets": data['datasets'][0],
                "revision_of": data['revision_of'][0]
            }
            if key == "links" or key == "datasets":
                val = ",".join(val)
            new_data[key] = val
            
        with tiledb.open(self.location + "/collections", "w") as A:
            A[id] = new_data

    def publish_collection(self, id):
        self.edit_collection(id, "published", True)

    def add_dataset(self, coll_id, name, artifact_id):
        id = new_id()
        datasets = self.get_attribute(coll_id, "datasets")
        datasets.append(id)
        self.edit_collection(coll_id, "datasets", datasets)
        with tiledb.open(self.location + "/datasets", "w") as A:
            A[id] = {
                "name": name,
                "artifact_id": artifact_id
            }
        return id

    def get_dataset(self, id):
        with tiledb.open(self.location + "/datasets", 'r') as A:
            return A[id]

    def edit_dataset(self, id, key, val):
        new_data = None
        with tiledb.open(self.location + "/datasets", "r") as A:
            data = A[id]
            new_data = {
                "name": data['name'][0],
                "artifact_id": data['artifact_id'][0]
            }
            new_data[key] = val
            
        with tiledb.open(self.location + "/datasets", "w") as A:
            A[id] = new_data

    def delete_dataset(self, coll_id: str, dataset_id: str):
        datasets = self.get_attribute(coll_id, "datasets")
        datasets.remove(dataset_id)
        self.edit_collection(coll_id, "datasets", datasets)
    
    def get_datasets(self, coll_id):
        ids = self.get_attribute(coll_id, "datasets")
        data = []
        with tiledb.open(self.location + "/datasets", 'r') as A:
            for id in ids:
                data.append(A[id])
        return data

    # TODO: maybe we also want the ability to read and revert based on timestamp, not just number of versions
    def read_collection_history(self, id, steps_back):
        fragments_info = tiledb.array_fragments(self.location + "/collections")
        if steps_back > len(fragments_info):
            raise IndexError("too many steps back in time")
        steps_idx = len(fragments_info) - steps_back - 1
        times = fragments_info.timestamp_range[steps_idx]

        with tiledb.open(self.location + "/collections", 'r', timestamp=times) as A:
            return A[id]

    def revert_collection_history(self, id, steps_back):
        fragments_info = tiledb.array_fragments(self.location + "/collections")
        if steps_back > len(fragments_info):
            raise IndexError("too many steps back in time")
        steps_idx = len(fragments_info) - steps_back - 1
        times = fragments_info.timestamp_range[steps_idx]

        data = None
        with tiledb.open(self.location + "/collections", 'r', timestamp=times) as A:
            data = A[id]

        # overwrite with old data
        new_data = {
            "published": data['published'][0],
            "name": data['name'][0],
            "description": data['description'][0],
            "contact_name": data['contact_name'][0],
            "contact_email": data['contact_email'][0],
            "links": data['links'][0], # TODO: figure out actual list support in TileDB
            "owner": data['owner'][0],
            "datasets": data['datasets'][0],
            "revision_of": data['revision_of'][0]
        }
        with tiledb.open(self.location + "/collections", 'w') as A:
            A[id] = new_data

    def read_dataset_history(self, id, steps_back):
        fragments_info = tiledb.array_fragments(self.location + "/datasets")
        if steps_back > len(fragments_info):
            raise IndexError("too many steps back in time")
        steps_idx = len(fragments_info) - steps_back - 1
        times = fragments_info.timestamp_range[steps_idx]

        with tiledb.open(self.location + "/datasets", 'r', timestamp=times) as A:
            return A[id]

    def revert_dataset_history(self, id, steps_back):
        fragments_info = tiledb.array_fragments(self.location + "/datasets")
        if steps_back > len(fragments_info):
            raise IndexError("too many steps back in time")
        steps_idx = len(fragments_info) - steps_back - 1
        times = fragments_info.timestamp_range[steps_idx]

        data = None
        with tiledb.open(self.location + "/datasets", 'r', timestamp=times) as A:
            data = A[id]

        # overwrite with old data
        new_data = {
            "name": data['name'][0],
            "artifact_id": data['artifact_id'][0]
        }
        with tiledb.open(self.location + "/datasets", 'w') as A:
            A[id] = new_data

    def create_revision(self, coll_id):
        id = new_id()
        coll = self.get_collection(coll_id)
        with tiledb.open(self.location + "/collections", "w") as A:
            A[id] = {
                "published": False,
                "name": coll['name'][0],
                "description": coll['description'][0],
                "contact_name": coll['contact_name'][0],
                "contact_email": coll['contact_email'][0],
                "links": coll['links'][0], 
                "owner": coll['owner'][0],
                "datasets": coll['datasets'][0],
                "revision_of": coll_id
            }
        return id

    def publish_revision(self, id):
        # get revision data
        revision = self.get_collection(id)
        new_data = {
            "published": False,
            "name": revision['name'][0],
            "description": revision['description'][0],
            "contact_name": revision['contact_name'][0],
            "contact_email": revision['contact_email'][0],
            "links": revision['links'][0],
            "owner": revision['owner'][0],
            "datasets": revision['datasets'][0],
            "revision_of": ""
        }
        # write data to existing revision_of collection
        revision_of = revision["revision_of"][0]
        with tiledb.open(self.location + "/collections", "w") as A:
            A[revision_of] = new_data
        # delete the revision
        self.delete_collection(id)

    def delete_collection(self, id):
        self.edit_collection(id, "published", -1)

    # TODO: functions for handling dataset artifacts (non-TileDB files)


    # TODO: functions to fit the APIs' requirements
