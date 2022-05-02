#!/usr/bin/env python3
import logging
import subprocess

from backend.corpora.dataset_processing.process import (
    get_bucket_prefix,
)

import time

import numpy as np
import psutil
import tiledb


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def process(dataset_id: str, cellxgene_bucket: str, prefix=None, dry_run=True):
    """
    Converts an existing CXG to the new, performance optimized schema
    :param dataset_id: The id of the dataset to reprocess.
    :param cellxgene_bucket: Name of the S3 bucket where the artifact is.
    :param prefix: If specified, will
    :param dry_run:
    :return:
    """

    # Download the cxg from the bucket
    if prefix is not None:
        object_key = f"{prefix}{dataset_id}"
    else:
        dataset_path = get_bucket_prefix(dataset_id)
        object_key = f"{dataset_path}.cxg"
    path = f"s3://{cellxgene_bucket}/{object_key}/X"

    logger.info(f"Processing dataset at path {path}, dry run {dry_run}")

    local_path = "/cxg"
    logger.info(f"Downloading {path} to {local_path}/X_old")

    download_command = ["aws", "s3", "sync", path, f"{local_path}/X_old"]
    # Let errors fail the pipeline
    subprocess.run(download_command, check=True)

    params = {
        "kind": "auto",
        "obs_extent": 256,
        "var_extent": 2048,
        "cell_order": "row",
        "tile_order": "col",
        "capacity": 128000,
        "compression": 22,
        "target_array": "X_new",
        "source_array": "X_old",
        "sparse_threshold": 25.0,
    }

    compute(cxg=local_path, **params)

    logger.info(f"Dataset at {path} computed successfully")

    if not dry_run:
        upload_command = ["aws", "s3", "sync", "--delete", f"{local_path}/X_new", path]
        subprocess.run(upload_command, check=True)

    # Cleanup
    logger.info("Cleaning up local files")
    import shutil

    shutil.rmtree(f"{local_path}/X_old")
    shutil.rmtree(f"{local_path}/X_new")


def compute(**kwargs):
    """
    Computes the evolved cxg from source_array and saves it to the target_array destination
    :param cxg:
    :param obs_extent:
    :param var_extent:
    :param source_array:
    :param target_array:
    :param kind:
    :param sparse_threshold:
    :param compression:
    :param cell_order:
    :param tile_order:
    :param capacity:
    """

    cxg = kwargs["cxg"]
    obs_extent = kwargs["obs_extent"]
    var_extent = kwargs["var_extent"]
    source_array = kwargs["source_array"]
    target_array = kwargs["target_array"]
    kind = kwargs["kind"]
    sparse_threshold = kwargs["sparse_threshold"]
    compression = kwargs["compression"]
    cell_order = kwargs["cell_order"]
    tile_order = kwargs["tile_order"]
    capacity = kwargs["capacity"]

    X_extent = [obs_extent, var_extent]

    ctx = create_fast_ctx(
        {
            "py.init_buffer_bytes": 4 * 1024**3,
            "sm.tile_cache_size": 1 * 1024**3,
        }
    )

    if kind == "auto":
        """
        Determine target encoding (sparse or dense) from the data.
        The "magic" sparsity threshold should in sync with
        backend/corpora/dataset_processing/process.py::make_cxg()
        """
        if not 0 < sparse_threshold <= 100:
            logger.info("Sparse threhold must be in range (0, 100]")
            return
        kind = choose_X_encoding(cxg, source_array, sparse_threshold, ctx)

    out_sparse = kind == "sparse"

    st = time.time()
    with tiledb.scope_ctx(ctx):
        with tiledb.open(f"{cxg}/{source_array}", "r") as old_X:
            create_new_X(
                cxg,
                target_array,
                kind,
                compression,
                X_extent,
                old_X.schema,
                cell_order,
                tile_order,
                capacity,
            )
            logger.info("created, starting to read...")
            in_sparse = old_X.schema.sparse
            with tiledb.open(f"{cxg}/{target_array}", "w") as new_X:
                logger.info(f"sparsity analysis: {in_sparse} {out_sparse} {old_X.schema.sparse} {new_X.schema.sparse}")
                if in_sparse and out_sparse:
                    logger.info("sparse->sparse")
                    i, chunk = 0, 200_000
                    while i < old_X.shape[0]:
                        logger.info(f"Chunk {i}")
                        dat = old_X[i : i + chunk]
                        logger.info("Got dat")
                        new_X[dat["obs"], dat["var"]] = dat[""]
                        i += chunk
                elif in_sparse and not out_sparse:
                    raise ValueError("sparse->dense not supported")
                elif not in_sparse and out_sparse:
                    logger.info("dense->sparse")
                    i, chunk = 0, 20_000
                    while i < old_X.shape[0]:
                        logger.info(f"Chunk {i}")
                        dat = old_X[i : i + chunk]
                        datnz = np.nonzero(dat)
                        new_X[datnz[0] + i, datnz[1]] = dat[datnz[0], datnz[1]]
                        i += chunk
                else:
                    logger.info("dense->dense")
                    new_X[:] = old_X[:]

            logger.info("consolidating...")
            tiledb.consolidate(f"{cxg}/{target_array}")
            tiledb.vacuum(f"{cxg}/{target_array}")

        logger.info(f"Took {time.time() - st} s")

    logger.info("done.")


def choose_X_encoding(cxg, source_array, sparse_threshold, ctx):
    logger.info("Reading array to determine sparsity...")
    with tiledb.scope_ctx(ctx):
        with tiledb.open(f"{cxg}/{source_array}", "r") as X:
            nnz = 0
            i, chunk = 0, 200_000
            while i < X.shape[0]:
                logger.info(f"Chunk {i}")
                dat = X.query(dims=[])[i : i + chunk]

                if X.schema.sparse:
                    nnz += len(dat[""])
                else:
                    nnz += np.count_nonzero(dat)
                i += chunk

            size = X.shape[0] * X.shape[1]
            logger.info(f"sparsity={1.0 - nnz / size}, nnz={nnz}, extent=({X.shape[0]}, {X.shape[1]})")
            if 1.0 - nnz / size >= sparse_threshold / 100.0:
                return "sparse"

            return "dense"


def create_new_X(cxg, target_array, kind, compression, X_extent, old_schema, cell_order, tile_order, capacity):
    is_sparse = kind == "sparse"
    logger.info(f"create_new_X: {is_sparse}")
    attr_filters = tiledb.FilterList([tiledb.ZstdFilter(level=compression)])
    dim_filters = tiledb.FilterList([tiledb.ByteShuffleFilter(), tiledb.ZstdFilter(level=compression)])
    old_dims = [old_schema.domain.dim(d) for d in range(old_schema.domain.ndim)]
    old_attr = old_schema.attr(0)

    tiledb.Array.create(
        f"{cxg}/{target_array}",
        tiledb.ArraySchema(
            domain=tiledb.Domain(
                [
                    tiledb.Dim(
                        name="obs",
                        domain=old_dims[0].domain,
                        tile=min(old_dims[0].domain[1], X_extent[0]),
                        dtype=old_dims[0].dtype,
                        filters=dim_filters,
                    ),
                    tiledb.Dim(
                        name="var",
                        domain=old_dims[1].domain,
                        tile=min(old_dims[1].domain[1], X_extent[1]),
                        dtype=old_dims[1].dtype,
                        filters=dim_filters,
                    ),
                ]
            ),
            sparse=is_sparse,
            allows_duplicates=True if is_sparse else False,
            attrs=[tiledb.Attr(name="", dtype=old_attr.dtype, filters=attr_filters)],
            cell_order=f"{cell_order}-major",
            tile_order=f"{tile_order}-major",
            capacity=capacity if is_sparse else 0,
        ),
    )


def create_ctx(config: dict = {}) -> tiledb.Ctx:
    cfg = tiledb.Config(config)
    ctx = tiledb.Ctx(config=cfg)
    return ctx


def frac_mem(f):
    mem_size = psutil.virtual_memory().total
    return int(f * mem_size) // (1024**2) * (1024**2)


def fast_config(config_overrides: dict = {}) -> dict:
    # consolidation buffer heuristic to prevent thrashing: total_mem/io_concurrency_level, rounded to GB
    io_concurrency_level = int(tiledb.Config()["sm.io_concurrency_level"])
    consolidation_buffer_size = (
        (int(frac_mem(0.1) / io_concurrency_level) + (1024**3 - 1)) // (1024**3) * (1024**3)
    )
    config = {
        "py.init_buffer_bytes": 16 * 1024**3,
        "sm.tile_cache_size": frac_mem(0.5),
        "sm.consolidation.buffer_size": consolidation_buffer_size,
    }
    config.update(config_overrides)
    return config


def create_fast_ctx(config_overrides: dict = {}) -> tiledb.Ctx:
    return create_ctx(fast_config(config_overrides))
