import os
import shutil
import numpy as np
import tiledb
import uuid
import ast
import time

"""
SCHEMA WITH TILEDB VERSIONING:

s3://single-cell-corpus/

    metadata/ [GROUP]
        collections [ARRAY]
            - uuid
            - owner
            - visibility
            - revision_of
            - datasets
            - other metadata
        datasets [ARRAY]
            - uuid
            - name
            - other metadata

    datasets/dataset_id/
        artifacts/
        uuid.h5ad
        uuid.rds
        uuid.cxg
    wmg/<SNAPSHOT>/
    soma/<SNAPSHOT>/
"""


class Utils:
    # TODO: config this somewhere, this is for testing only
    location = "tests/unit/backend/fixtures/test_tiledb/metadata"

    attrs = {
        "collections": [
            # API
            "name",
            "description",
            "contact_name",
            "contact_email",
            "links",
            "datasets",
            "curator_name",
            "created_at",
            "updated_at",
            "publisher_metadata",
            # internal TileDB schema
            "visibility",
            "revision_of",
            "owner",
        ],
        "datasets": [
            "x_approximate_distribution",
            "x_normalization",
            "assay",
            "cell_count",
            "cell_type",
            "development_stage",
            "disease",
            "ethnicity",
            "is_primary_data",
            "name",
            "organism",
            "sex",
            "tissue",
            "explorer_url",
            "processing_status",
            "dataset_assets"
        ]
    }
    attrs_to_parse = {
        "collections":   [
            "links",
            "datasets",
            "publisher_metadata"
        ],
        "datasets": [
            "assay",
            "cell_type",
            "development_stage",
            "disease",
            "ethnicity",
            "organism",
            "sex",
            "tissue",
            "dataset_assets",
            "processing_status"
        ]
    }

    empty_dataset = {
        "x_approximate_distribution": "PROCESSING",
        "x_normalization": "",
        "assay": [],
        "cell_count": 0,
        "cell_type": [],
        "development_stage": [],
        "disease": [],
        "ethnicity": [],
        "dataset_assets": [],
        "is_primary_data": "PROCESSING",
        "name": "",
        "organism": [],
        "sex": [],
        "tissue": [],
        "explorer_url": "",
        "processing_status": {}
    }

    empty_collection = {
        "name": "",
        "description": "",
        "owner": "",
        "contact_name": "",
        "contact_email": "",
        "curator_name": "",
        "links": [],
        "publisher_metadata": {},
        "datasets": [],
        "created_at": 0,
        "updated_at": 0,
        "visibility": "PRIVATE",
        "revision_of": ""
    }

    @staticmethod
    def new_id():
        return uuid.uuid4().hex

    # TODO: figure out actual list and dict support in TileDB
    @staticmethod
    def parse_stored_data(data: dict, array: str) -> dict:
        for a in Utils.attrs[array]:
            if a in data:
                t = type(data[a])
                if t == np.ndarray:
                    if len(data[a]) > 0:
                        data[a] = data[a][0]  # for some reason some fields get stored as arrays in TileDB
                    else:
                        if array == "collections":
                            data[a] = Utils.empty_collection[a]
                        else:
                            data[a] = Utils.empty_dataset[a]
                        if a in Utils.attrs_to_parse[array]:
                            data[a] = str(data[a])
                if a in Utils.attrs_to_parse[array]:
                    data[a] = ast.literal_eval(data[a])
                t = type(data[a])
                if t == np.float32:
                    data[a] = float(data[a])
                elif t == np.int32:
                    data[a] = int(data[a])
        if type(data['id']) == np.ndarray and len(data['id']) > 0:
            data['id'] = data['id'][0].decode("utf-8") # TileDB stores the id index as byte string
        elif len(data['id']) > 0:
            data['id'] = data['id'].decode("utf-8")
        return data

    @staticmethod
    def pack_input_data(data: dict, array: str) -> dict:
        for a in Utils.attrs_to_parse[array]:
            if a in data:
                data[a] = str(data[a])
        return data


class TileDBData:
    @staticmethod
    def init_db(location=Utils.location):
        """Create a local TileDB group and arrays according to our schema."""
        # create group
        if os.path.exists(location):
            shutil.rmtree(location)

        tiledb.group_create(location)

        # create collections array
        # TODO: figure out ideal domain, tile, and number of dimensions
        dim1 = tiledb.Dim(name="id", domain=(None, None), tile=2, dtype="S0")
        dom = tiledb.Domain(dim1)

        a1 = tiledb.Attr(name="owner", dtype="U1")
        a2 = tiledb.Attr(name="visibility", dtype="U1")  # DELETED, PRIVATE, PUBLIC
        a3 = tiledb.Attr(name="name", dtype="U1")
        a4 = tiledb.Attr(name="description", dtype="U1")
        a5 = tiledb.Attr(name="contact_name", dtype="U1")
        a6 = tiledb.Attr(name="contact_email", dtype="U1")
        a7 = tiledb.Attr(name="links", dtype="U1")
        a8 = tiledb.Attr(name="datasets", dtype="U1")
        a9 = tiledb.Attr(name="revision_of", dtype="U1")
        a10 = tiledb.Attr(name="curator_name", dtype="U1")
        a11 = tiledb.Attr(name="created_at", dtype=np.float32)
        a12 = tiledb.Attr(name="updated_at", dtype=np.float32)
        a13 = tiledb.Attr(name="publisher_metadata", dtype="U1")
        attrs = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13]

        schema = tiledb.ArraySchema(domain=dom, sparse=True, attrs=attrs)
        array = location + "/collections"
        tiledb.Array.create(array, schema)

        # create datasets array
        a1 = tiledb.Attr(name="x_approximate_distribution", dtype="U1")
        a2 = tiledb.Attr(name="x_normalization", dtype="U1")
        a3 = tiledb.Attr(name="cell_count", dtype=np.int32)
        a4 = tiledb.Attr(name="cell_type", dtype="U1")
        a5 = tiledb.Attr(name="development_stage", dtype="U1")
        a6 = tiledb.Attr(name="disease", dtype="U1")
        a7 = tiledb.Attr(name="ethnicity", dtype="U1")
        a8 = tiledb.Attr(name="is_primary_data", dtype="U1")
        a9 = tiledb.Attr(name="name", dtype="U1")
        a10 = tiledb.Attr(name="organism", dtype="U1")
        a11 = tiledb.Attr(name="sex", dtype="U1")
        a12 = tiledb.Attr(name="tissue", dtype="U1")
        a13 = tiledb.Attr(name="explorer_url", dtype="U1")
        a14 = tiledb.Attr(name="processing_status", dtype="U1")
        a15 = tiledb.Attr(name="assay", dtype="U1")
        a16 = tiledb.Attr(name="dataset_assets", dtype="U1")
        attrs = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16]

        schema = tiledb.ArraySchema(domain=dom, sparse=True, attrs=attrs)
        array = location + "/datasets"
        tiledb.Array.create(array, schema)

    @staticmethod
    def destroy_db(location):
        """Delete our local TileDB group."""
        if os.path.exists(location):
            shutil.rmtree(location)

    def __init__(self, location=Utils.location):
        self.location = location

    def create_collection(self, metadata: dict = {}):
        """Creates a collection using provided data."""
        id = Utils.new_id()
        return self._create_collection_custom_id(id, metadata)

    def _create_collection_custom_id(self, id: str = "", metadata: dict = {}):
        """Internal function to create a collection with a custom id. """
        data = Utils.empty_collection
        for a in Utils.attrs['collections']:
            if a in metadata and metadata[a]:
                data[a] = metadata[a]
        data["created_at"] = time.time()
        data["updated_at"] = time.time()
        with tiledb.open(self.location + "/collections", "w") as A:
            A[id] = Utils.pack_input_data(data, "collections")
        return id

    def get_collection(self, id):
        """Gets a collection by its id"""
        with tiledb.open(self.location + "/collections", 'r') as A:
            res = Utils.parse_stored_data(dict(A[id]), "collections")
            n = res['id']
            return res if (type(n) == str and len(n) > 0) or (type(n) == np.ndarray and n.size > 0) else None

    def get_published_collections(self):
        """Get all public collections"""
        with tiledb.open(self.location + "/collections", mode="r") as A:
            # TODO: try query conditions for efficiency
            df = A.df[:]
            df = df[(df['visibility'] == "PUBLIC")]
            res = df.to_dict("records")
            for i in range(len(res)):
                res[i] = Utils.parse_stored_data(res[i], "collections")
            return res

    def get_all_collections(self):
        """Get all collections"""
        with tiledb.open(self.location + "/collections", mode="r") as A:
            df = A.df[:]
            res = df.to_dict("records")
            for i in range(len(res)):
                res[i] = Utils.parse_stored_data(res[i], "collections")
            return res

    def get_published_datasets(self):
        """Get all datasets belonging to a public collection"""
        colls = self.get_published_collections()
        dataset_ids = []
        for coll in colls:
            d = coll['datasets']
            dataset_ids.extend(d if d else [])
        with tiledb.open(self.location + "/datasets", mode="r") as A:
            df = A.df[:]
            df = df[(df['id'].isin(dataset_ids))]
            res = df.to_dict("records")
            for i in range(len(res)):
                res[i] = Utils.parse_stored_data(res[i], "datasets")
            return res

    def get_attribute(self, id, attr):
        """Get the data stored in one field of a specific collection"""
        coll = self.get_collection(id)
        data = coll[attr]
        return data

    def edit_collection(self, id, key, val):
        """Update the data stored in one field of a specific collection"""
        new_data = None
        with tiledb.open(self.location + "/collections", "r") as A:
            data = A[id]
            new_data = {}
            for attr in Utils.attrs["collections"]:
                new_data[attr] = data[attr][0] if len(data[attr]) > 0 else Utils.empty_collection[attr]
            new_data[key] = val
            new_data["updated_at"] = time.time()
            new_data = Utils.pack_input_data(new_data, "collections")

        with tiledb.open(self.location + "/collections", "w") as A:
            A[id] = new_data

    def publish_collection(self, id: str):
        """Set a collection's visibility to public"""
        self.edit_collection(id, "visibility", "PUBLIC")

    def add_dataset(self, coll_id: str, metadata: dict):
        """Add a dataset to a collection and to the datasets array using the data from the user's shared URL"""
        id = Utils.new_id()
        return self._add_dataset_custom_id(id, coll_id, metadata)

    def _add_dataset_custom_id(self, id: str, coll_id: str, metadata: dict):
        """Add a dataset with a custom id."""
        datasets = self.get_attribute(coll_id, "datasets")
        if id not in datasets:
            datasets.append(id)
        self.edit_collection(coll_id, "datasets", datasets)

        data = Utils.empty_dataset
        for a in Utils.attrs['datasets']:
            if a in metadata and metadata[a]:
                data[a] = metadata[a]

        with tiledb.open(self.location + "/datasets", "w") as A:
            A[id] = Utils.pack_input_data(data, "datasets")
        return id

    def get_dataset(self, id: str):
        """Get a dataset by its id"""
        with tiledb.open(self.location + "/datasets", 'r') as A:
            res = Utils.parse_stored_data(dict(A[id]), "datasets")
            n = res['id']
            return res if (type(n) == str and len(n) > 0) or (type(n) == np.ndarray and n.size > 0) else None

    def edit_dataset(self, id, key, val):
        """Update the data in one field of a specific dataset"""
        new_data = None
        with tiledb.open(self.location + "/datasets", "r") as A:
            data = A[id]
            new_data = {}
            for attr in Utils.attrs["datasets"]:
                new_data[attr] = data[attr][0] if len(data[attr]) > 0 else Utils.empty_dataset[attr]
            new_data[key] = val
            new_data = Utils.pack_input_data(new_data, "datasets")

        with tiledb.open(self.location + "/datasets", "w") as A:
            A[id] = new_data

    def bulk_edit_dataset(self, id: str, metadata: str):
        """Replace all data in a dataset by its id"""
        with tiledb.open(self.location + "/datasets", "w") as A:
            A[id] = Utils.pack_input_data(metadata, "datasets")

    def delete_dataset(self, coll_id: str, dataset_id: str):
        """Remove a dataset from a collection"""
        datasets = self.get_attribute(coll_id, "datasets")
        datasets.remove(dataset_id)
        self.edit_collection(coll_id, "datasets", datasets)

    def replace_dataset(self, dataset_id: str, assets: list):
        self.edit_dataset(dataset_id, "dataset_assets", assets)

    def get_datasets(self, coll_id):
        """Return all the datasets belonging to a specific collection"""
        ids = self.get_attribute(coll_id, "datasets")
        data = []
        with tiledb.open(self.location + "/datasets", 'r') as A:
            for id in ids:
                data.append(Utils.parse_stored_data(dict(A[id]), "datasets"))
        return data

    # TODO: maybe we also want the ability to read and revert based on timestamp, not just number of versions
    def read_collection_history(self, id, steps_back):
        """Get a collection by its id some specific number of writes ago"""
        fragments_info = tiledb.array_fragments(self.location + "/collections")
        if steps_back > len(fragments_info):
            raise IndexError("too many steps back in time")
        steps_idx = len(fragments_info) - steps_back - 1
        times = fragments_info.timestamp_range[steps_idx]

        with tiledb.open(self.location + "/collections", 'r', timestamp=times) as A:
            return Utils.parse_stored_data(dict(A[id]), "collections")

    def revert_collection_history(self, id, steps_back):
        """Revert a collection by its id to the state it was in a specific number of writes ago"""
        fragments_info = tiledb.array_fragments(self.location + "/collections")
        if steps_back > len(fragments_info):
            raise IndexError("too many steps back in time")
        steps_idx = len(fragments_info) - steps_back - 1
        times = fragments_info.timestamp_range[steps_idx]

        data = None
        with tiledb.open(self.location + "/collections", 'r', timestamp=times) as A:
            data = A[id]

        # overwrite with old data
        new_data = {}
        for attr in Utils.attrs["collections"]:
            new_data[attr] = data[attr][0]
        new_data["updated_at"] = time.time()
        new_data = Utils.pack_input_data(new_data, "collections")
        with tiledb.open(self.location + "/collections", 'w') as A:
            A[id] = new_data

    def read_dataset_history(self, id, steps_back):
        """Get a dataset by its id some specific number of writes ago"""
        fragments_info = tiledb.array_fragments(self.location + "/datasets")
        if steps_back > len(fragments_info):
            raise IndexError("too many steps back in time")
        steps_idx = len(fragments_info) - steps_back - 1
        times = fragments_info.timestamp_range[steps_idx]

        with tiledb.open(self.location + "/datasets", 'r', timestamp=times) as A:
            return Utils.parse_stored_data(dict(A[id]), "datasets")

    def revert_dataset_history(self, id, steps_back):
        """Revert a dataset by its id to the state it was in a specific number of writes ago"""
        fragments_info = tiledb.array_fragments(self.location + "/datasets")
        if steps_back > len(fragments_info):
            raise IndexError("too many steps back in time")
        steps_idx = len(fragments_info) - steps_back - 1
        times = fragments_info.timestamp_range[steps_idx]

        data = None
        with tiledb.open(self.location + "/datasets", 'r', timestamp=times) as A:
            data = A[id]

        # overwrite with old data
        new_data = {}
        for a in Utils.attrs['datasets']:
            new_data[a] = data[a][0]

        with tiledb.open(self.location + "/datasets", 'w') as A:
            A[id] = new_data

    def create_revision(self, coll_id: str):
        """Start a revision of an existing collection by its id"""
        id = Utils.new_id()
        return self._create_revision_custom_id(id, coll_id)

    def _create_revision_custom_id(self, id: str, coll_id: str):
        """Start a revision with a custom id"""
        coll = self.get_collection(coll_id)
        data = {}
        for attr in Utils.attrs["collections"]:
            data[attr] = coll[attr]
        data['revision_of'] = coll_id
        data["visibility"] = "PRIVATE"
        with tiledb.open(self.location + "/collections", "w") as A:
            A[id] = Utils.pack_input_data(data, "collections")
        return id

    def publish_revision(self, id):
        """Publish a revision by its id"""
        # get revision data
        revision = self.get_collection(id)
        data = {}
        for attr in Utils.attrs["collections"]:
            data[attr] = revision[attr]
        data['visibility'] = "PUBLIC"
        data['revision_of'] = ""
        data['updated_at'] = time.time()
        # write data to existing revision_of collection
        revision_of = revision["revision_of"]
        with tiledb.open(self.location + "/collections", "w") as A:
            A[revision_of] = Utils.pack_input_data(data, "collections")
        # delete the revision
        self.delete_collection(id)

    def delete_collection(self, id):
        """Mark a collection as deleted by its id"""
        self.edit_collection(id, "visibility", "DELETED")

