import logging

import tiledb
from scipy import sparse

logger = logging.getLogger(__name__)


class TileDBSparseArrayWriteWrapper:
    def __init__(self, uri, *, ctx=None):
        self.uri = uri
        self.ctx = {}
        self.ctx.update(**ctx or {})

    def __setitem__(self, k: tuple[slice, ...], v: sparse.spmatrix):
        with tiledb.scope_ctx(self.ctx):
            row_slice, col_slice = k
            row_offset = row_slice.start if row_slice.start is not None else 0
            col_offset = col_slice.start if col_slice.start is not None else 0
            v_coo = v.tocoo()
            tiledb_array = tiledb.open(self.uri, mode="w")
            tiledb_array[v_coo.row + row_offset, v_coo.col + col_offset] = v.data
