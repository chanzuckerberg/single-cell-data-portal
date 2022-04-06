import argparse
import time

import numpy as np
import psutil
import tiledb

"""
Re-master an existing "X" tiledb array in a CXG
Creates a new array with the same data but a different tdb schema
"""

X_extent = None
X_name = None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cxg", type=str, help="cxg path")
    parser.add_argument(
        "--kind", choices=["auto", "sparse", "dense"], type=str, help="target: auto, sparse or dense", default="auto"
    )
    parser.add_argument("--obs", type=int, help="obs extent", default=256)
    parser.add_argument("--var", type=int, help="var extent", default=2048)
    parser.add_argument("--cell-order", choices=["row", "col"], type=str, default="row")
    parser.add_argument("--tile-order", choices=["row", "col"], type=str, default="col")
    parser.add_argument("--capacity", type=int, default=128000)
    parser.add_argument("--compression", type=int, default=22)
    parser.add_argument("--name", type=str, default="X_new")
    parser.add_argument("--arr", type=str, default="X")
    parser.add_argument("--sparse-threshold", type=float, default=25.0, help="sparse encoding threshold")
    args = parser.parse_args()

    global X_extent, X_name
    X_extent = [args.obs, args.var]
    X_name = args.name

    ctx = create_fast_ctx(
        {
            "py.init_buffer_bytes": 32 * 1024 ** 3,
            "sm.tile_cache_size": 2 * 1024 ** 3,
        }
    )

    if args.kind == "auto":
        """
        Determine target encoding (sparse or dense) from the data.
        The "magic" sparsity threshold should in sync with
        backend/corpora/dataset_processing/process.py::make_cxg()
        """
        if not 0 < args.sparse_threshold <= 100:
            print("Sparse threhold must be in range (0, 100]")
            return
        args.kind = choose_X_encoding(args, ctx)

    out_sparse = args.kind == "sparse"

    st = time.time()
    with tiledb.scope_ctx(ctx):
        with tiledb.open(f"{args.cxg}/{args.arr}", "r") as old_X:
            create_new_X(args, old_X.schema)
            print("created, starting to read...")
            in_sparse = old_X.schema.sparse
            with tiledb.open(f"{args.cxg}/{X_name}", "w") as new_X:
                if in_sparse and out_sparse:
                    print("sparse->sparse")
                    i, chunk = 0, 200_000
                    while i < old_X.shape[0]:
                        print(f"Chunk {i}")
                        dat = old_X[i : i + chunk]
                        print(f"Got dat")
                        new_X[dat["obs"], dat["var"]] = dat[""]
                        i += chunk
                elif in_sparse and not out_sparse:
                    raise ValueError("sparse->dense not supported")
                elif not in_sparse and out_sparse:
                    print("dense->sparse")
                    i, chunk = 0, 20_000
                    while i < old_X.shape[0]:
                        print(f"Chunk {i}")
                        dat = old_X[i : i + chunk]
                        datnz = np.nonzero(dat)
                        new_X[datnz[0] + i, datnz[1]] = dat[datnz[0], datnz[1]]
                        i += chunk
                else:
                    print("dense->dense")
                    new_X[:] = old_X[:]

            print("consolidating...")
            tiledb.consolidate(f"{args.cxg}/{X_name}")
            tiledb.vacuum(f"{args.cxg}/{X_name}")

        print(f"Took {time.time() - st} s")

    print("done.")


def choose_X_encoding(args, ctx):
    print("Reading array to determine sparsity...")
    with tiledb.scope_ctx(ctx):
        with tiledb.open(f"{args.cxg}/{args.arr}", "r") as X:
            nnz = 0
            i, chunk = 0, 200_000
            while i < X.shape[0]:
                print(f"Chunk {i}")
                dat = X.query(dims=[])[i : i + chunk]

                if X.schema.sparse:
                    nnz += len(dat[""])
                else:
                    nnz += np.count_nonzero(dat)
                i += chunk

            size = X.shape[0] * X.shape[1]
            print(f"sparsity={1.0 - nnz / size}, nnz={nnz}, extent=({X.shape[0]}, {X.shape[1]})")
            if 1.0 - nnz / size >= args.sparse_threshold / 100.0:
                return "sparse"

            return "dense"


def create_new_X(args, old_schema):
    is_sparse = args.kind == "sparse"
    attr_filters = tiledb.FilterList([tiledb.ZstdFilter(level=args.compression)])
    dim_filters = tiledb.FilterList([tiledb.ByteShuffleFilter(), tiledb.ZstdFilter(level=args.compression)])
    old_dims = [old_schema.domain.dim(d) for d in range(old_schema.domain.ndim)]
    old_attr = old_schema.attr(0)

    tiledb.Array.create(
        f"{args.cxg}/{X_name}",
        tiledb.ArraySchema(
            domain=tiledb.Domain(
                [
                    tiledb.Dim(
                        name="obs",
                        domain=old_dims[0].domain,
                        tile=X_extent[0],
                        dtype=old_dims[0].dtype,
                        filters=dim_filters,
                    ),
                    tiledb.Dim(
                        name="var",
                        domain=old_dims[1].domain,
                        tile=X_extent[1],
                        dtype=old_dims[1].dtype,
                        filters=dim_filters,
                    ),
                ]
            ),
            sparse=is_sparse,
            allows_duplicates=True if is_sparse else False,
            attrs=[tiledb.Attr(name="", dtype=old_attr.dtype, filters=attr_filters)],
            cell_order=f"{args.cell_order}-major",
            tile_order=f"{args.tile_order}-major",
            capacity=args.capacity if is_sparse else 0,
        ),
    )


def create_ctx(config: dict = {}) -> tiledb.Ctx:
    cfg = tiledb.Config(config)
    ctx = tiledb.Ctx(config=cfg)
    return ctx


def frac_mem(f):
    mem_size = psutil.virtual_memory().total
    return int(f * mem_size) // (1024 ** 2) * (1024 ** 2)


def fast_config(config_overrides: dict = {}) -> dict:
    # consolidation buffer heuristic to prevent thrashing: total_mem/io_concurrency_level, rounded to GB
    io_concurrency_level = int(tiledb.Config()["sm.io_concurrency_level"])
    consolidation_buffer_size = (
        (int(frac_mem(0.1) / io_concurrency_level) + (1024 ** 3 - 1)) // (1024 ** 3) * (1024 ** 3)
    )
    config = {
        "py.init_buffer_bytes": 16 * 1024 ** 3,
        "sm.tile_cache_size": frac_mem(0.5),
        "sm.consolidation.buffer_size": consolidation_buffer_size,
    } | config_overrides
    config.update(config)
    return config


def create_fast_ctx(config_overrides: dict = {}) -> tiledb.Ctx:
    return create_ctx(fast_config(config_overrides))


if __name__ == "__main__":
    main()
