import json
import logging

import dask.array as da
import numpy as np
import pandas as pd
import tiledb
from cellxgene_schema.utils import get_matrix_format

from backend.common.constants import IS_SINGLE, UNS_SPATIAL_KEY
from backend.layers.processing.utils.atac import ATACDataProcessor
from backend.layers.processing.utils.dask_utils import TileDBSparseArrayWriteWrapper
from backend.layers.processing.utils.spatial import SpatialDataProcessor
from backend.layers.processing.utils.type_conversion_utils import get_dtype_and_schema_of_array


def convert_dictionary_to_cxg_group(cxg_container, metadata_dict, group_metadata_name="cxg_group_metadata", ctx=None):
    """
    Saves the contents of the dictionary to the CXG output directory specified.

    This function is primarily used to save metadata about a dataset to the CXG directory. At some point, tiledb will
    have support for metadata on groups at which point the utility of this function should be revisited. Until such
    feature exists, this function create an empty array and annotate that array.

    For more information, visit https://github.com/TileDB-Inc/TileDB-Py/issues/254.
    """

    array_name = f"{cxg_container}/{group_metadata_name}"

    # Because TileDB does not allow one to attach metadata directly to a CXG group, we need to have a workaround
    # where we create an empty array and attached the metadata onto to this empty array. Below we construct this empty
    # array.
    tiledb.from_numpy(array_name, np.zeros((1,)))

    with tiledb.open(array_name, mode="w", ctx=ctx) as metadata_array:
        for key, value in metadata_dict.items():
            metadata_array.meta[key] = value


def convert_uns_to_cxg_group(cxg_container, metadata_dict, dataset_version_id, group_metadata_name="uns", ctx=None):
    """
    Convert uns (unstructured) metadata to CXG output directory specified
    Generate deep zoom assets for spatial data
    """

    spatial_processor = SpatialDataProcessor()

    array_name = f"{cxg_container}/{group_metadata_name}"
    object_filtered = {}

    tiledb.from_numpy(array_name, np.zeros((1,)))

    with tiledb.open(array_name, mode="w", ctx=ctx) as metadata_array:
        for key, value in metadata_dict.items():
            if key == UNS_SPATIAL_KEY:
                for object_id, content in value.items():
                    if object_id != IS_SINGLE:
                        object_filtered = spatial_processor.filter_spatial_data(content, object_id)
                        spatial_processor.create_deep_zoom_assets(cxg_container, content, dataset_version_id)

                metadata_array.meta[key] = json.dumps(object_filtered)


def convert_coverage_to_cxg_array(cxg_container, metadata_dict, fragment_artifact_id, group_metadata_name, ctx):

    atac_processor = ATACDataProcessor(fragment_artifact_id, ctx)

    array_name = f"{cxg_container}/{group_metadata_name}"

    df_meta, cell_id_map = atac_processor.process_fragment_file(metadata_dict, array_name)

    print(df_meta)

    with tiledb.open(array_name, mode="w", ctx=ctx) as array:
        array.meta["cell_id_map"] = json.dumps(cell_id_map)
        cell_metadata_dict = df_meta.set_index("cell_name")["cell_type"].to_dict()
        array.meta["cell_metadata"] = json.dumps(cell_metadata_dict)

    # tiledb.consolidate(array_name, ctx=ctx)


def convert_dataframe_to_cxg_array(cxg_container, dataframe_name, dataframe, index_column_name, ctx):
    """
    Saves the contents of the dataframe to the CXG output directory specified.

    Current access patterns are oriented toward reading very large slices of the dataframe, one attribute at a time.
    Attribute data also tends to be (often) repetitive (bools, categories, strings). Given this, we use a large tile
    size (1000) and very aggressive compression levels.
    """

    tiledb_filter = tiledb.FilterList(
        [
            # Attempt aggressive compression as many of these dataframes are very repetitive strings, bools and
            # other non-float data.
            tiledb.ZstdFilter(level=22),
        ]
    )
    data = {}
    schema_hints = {}
    tdb_attrs = []

    for column_name, column_values in dataframe.items():
        # Cast 'in_tissue' column values as boolean to make it categorical
        # https://github.com/chanzuckerberg/single-cell-explorer/issues/841
        if column_name == "in_tissue":
            dtype, hints = get_dtype_and_schema_of_array(column_values.astype(bool))
        else:
            dtype, hints = get_dtype_and_schema_of_array(column_values)
        if "categories" in hints and len(hints.get("categories", [])) > 0.75 * dataframe.shape[0]:
            hints["type"] = "string"
            del hints["categories"]
            data[column_name] = column_values.to_numpy(dtype=dtype)
        elif "categories" in hints:
            cat = pd.Categorical(column_values.astype("str"))
            codes = cat.codes
            data[column_name] = codes
            categories = list(cat.categories)
            hints["categories"] = categories
            dtype = str(cat.codes.dtype)
        else:
            data[column_name] = column_values.to_numpy(dtype=dtype)

        tdb_attrs.append(tiledb.Attr(name=column_name, dtype=dtype, filters=tiledb_filter))
        schema_hints.update({column_name: hints})

    def create_dataframe_array(array_name, dataframe, attrs):
        domain = tiledb.Domain(
            tiledb.Dim(domain=(0, dataframe.shape[0] - 1), tile=min(dataframe.shape[0], 1000), dtype=np.uint32)
        )
        schema = tiledb.ArraySchema(
            domain=domain, sparse=False, attrs=attrs, cell_order="row-major", tile_order="row-major"
        )
        tiledb.DenseArray.create(array_name, schema)

    array_name = f"{cxg_container}/{dataframe_name}"

    create_dataframe_array(array_name, dataframe, tdb_attrs)

    with tiledb.open(array_name, mode="w", ctx=ctx) as array:
        schema_hints.update({"index": index_column_name})
        array[:] = data
        array.meta["cxg_schema"] = json.dumps(schema_hints)

    tiledb.consolidate(array_name, ctx=ctx)


def convert_ndarray_to_cxg_dense_array(ndarray_name, ndarray, ctx):
    """
    Saves contents of ndarray to the CXG output directory specified.

    Generally this function is used to convert dataset embeddings. Because embeddings are typically accessed with
    very large slices (or all of the embedding), they do not benefit from overly aggressive compression due to their
    format.  Given this, we use a large tile size (1000) but only default compression level.
    """

    def create_ndarray_array(ndarray_name, ndarray):
        filters = tiledb.FilterList([tiledb.ZstdFilter()])
        attrs = [tiledb.Attr(dtype=ndarray.dtype, filters=filters)]
        dimensions = [
            tiledb.Dim(
                domain=(0, ndarray.shape[dimension] - 1), tile=min(ndarray.shape[dimension], 1000), dtype=np.uint32
            )
            for dimension in range(ndarray.ndim)
        ]
        domain = tiledb.Domain(*dimensions)
        schema = tiledb.ArraySchema(
            domain=domain, sparse=False, attrs=attrs, capacity=1_000_000, cell_order="row-major", tile_order="row-major"
        )
        tiledb.DenseArray.create(ndarray_name, schema)

    create_ndarray_array(ndarray_name, ndarray)

    with tiledb.open(ndarray_name, mode="w", ctx=ctx) as array:
        array[:] = ndarray

    tiledb.consolidate(ndarray_name, ctx=ctx)


def convert_matrices_to_cxg_arrays(matrix_name: str, matrix: da.Array, encode_as_sparse_array: bool, ctx: tiledb.Ctx):
    """
    Converts a numpy array matrix into a TileDB SparseArray of DenseArray based on whether `encode_as_sparse_array`
    is true or not. Note that when the matrix is encoded as a SparseArray, it only writes the values that are
    nonzero. This means that if you count the number of elements in the SparseArray, it will not equal the total
    number of elements in the matrix, only the number of nonzero elements.
    """
    number_of_rows = matrix.shape[0]
    number_of_columns = matrix.shape[1]
    compression = 22

    logging.info(f"create {matrix_name}")
    dim_filters = tiledb.FilterList([tiledb.ByteShuffleFilter(), tiledb.ZstdFilter(level=compression)])
    attrs = [tiledb.Attr(dtype=np.float32, filters=tiledb.FilterList([tiledb.ZstdFilter(level=compression)]))]

    tiledb_obs_dim = tiledb.Dim(
        name="obs",
        domain=(0, number_of_rows - 1),
        tile=min(number_of_rows, 256),
        dtype=np.uint32,
        filters=dim_filters,
    )
    tiledb_var_dim = tiledb.Dim(
        name="var",
        domain=(0, number_of_columns - 1),
        tile=min(number_of_columns, 2048),
        dtype=np.uint32,
        filters=dim_filters,
    )
    domain = tiledb.Domain(tiledb_obs_dim, tiledb_var_dim)

    if encode_as_sparse_array:
        array_schema_params = dict(
            sparse=True,
            allows_duplicates=True,
            capacity=1024000,
        )
    else:
        array_schema_params = dict(
            sparse=False,
            allows_duplicates=False,
            capacity=0,
        )
    schema = tiledb.ArraySchema(
        domain=domain,
        attrs=attrs,
        cell_order="row-major",
        tile_order="col-major",
        **array_schema_params,
    )
    tiledb.Array.create(matrix_name, schema)

    if encode_as_sparse_array:
        matrix_write = TileDBSparseArrayWriteWrapper(matrix_name, ctx=ctx)
        matrix.store(matrix_write, lock=False, compute=True)
    else:
        # if matrix is a scipy sparse matrix but encode_as_sparse_array is False, convert to dense array
        if get_matrix_format(matrix) != "dense":
            matrix = matrix.map_blocks(
                lambda x: x.toarray().astype(np.float32), dtype=np.float32, meta=np.array([], dtype=np.float32)
            )
        elif matrix.dtype != np.float32:
            matrix = matrix.map_blocks(lambda x: x.astype(np.float32), dtype=np.float32)
        with tiledb.open(matrix_name, "w") as A:
            matrix.to_tiledb(A, storage_options={"ctx": ctx})
