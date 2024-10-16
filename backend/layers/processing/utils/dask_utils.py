import logging

import dask.distributed as dd
import numpy as np
import tiledb
from scipy import sparse

logger = logging.getLogger(__name__)


class TileDBSparseArrayWriteWrapper:
    def __init__(self, uri, *, ctx=None):

        self.uri = uri
        self.ctx = ctx or {}

    def __setitem__(self, k: tuple[slice, ...], v: sparse.spmatrix):
        with tiledb.scope_ctx(
            {
                "sm.io_concurrency_level": 1,
                "sm.compute_concurrency_level": 1,
            }
        ):
            row_slice, col_slice = k
            row_offset = row_slice.start if row_slice.start is not None else 0
            col_offset = col_slice.start if col_slice.start is not None else 0
            v_coo = v.tocoo()
            tiledb_array = tiledb.open(self.uri, mode="w")
            tiledb_array[v_coo.row + row_offset, v_coo.col + col_offset] = v.data


def write_dask_array_as_tiledb(output_path, X, *, filters=None):
    filters = filters or tiledb.FilterList([tiledb.ZstdFilter()])
    # Create output array
    attrs = [tiledb.Attr(dtype=X.dtype, filters=filters)]
    domain = tiledb.Domain(
        tiledb.Dim(
            name="obs",
            domain=(0, X.shape[0] - 1),
            tile=min(X.shape[0], 256),
            dtype=np.uint32,
            filters=filters,
        ),
        tiledb.Dim(
            name="var",
            domain=(0, X.shape[1] - 1),
            tile=min(X.shape[1], 2048),
            dtype=np.uint32,
            filters=filters,
        ),
    )
    schema = tiledb.ArraySchema(
        domain=domain,
        sparse=True,
        allows_duplicates=True,
        attrs=attrs,
        cell_order="row-major",
        tile_order="row-major",
        capacity=1024000,
    )
    tiledb.Array.create(output_path, schema)

    # Write
    X_write = TileDBSparseArrayWriteWrapper(output_path)
    X.store(X_write, lock=False, compute=True)


DASK_CLIENT = None


def start_dask_cluster():
    global DASK_CLIENT
    if DASK_CLIENT:
        return DASK_CLIENT
    logger.info("Starting Dask cluster")
    cluster = dd.LocalCluster()
    client = dd.Client(cluster)
    return client
