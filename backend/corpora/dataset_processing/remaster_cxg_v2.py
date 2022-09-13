#!/usr/bin/env python3
import logging
import subprocess
import json
import numpy as np
import pandas as pd
from backend.corpora.dataset_processing.common import get_bucket_prefix
import time

import psutil
import tiledb


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def process(dataset_id: str, cellxgene_bucket: str, prefix=None, dry_run=True):
    """
    Converts an existing CXG to the new schema with row- and column-oriented arrays
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
    obs_path = f"s3://{cellxgene_bucket}/{object_key}/obs"

    logger.info(f"Processing dataset at path {path}, dry run {dry_run}")

    local_path = "/cxg"
    logger.info(f"Downloading {path} to {local_path}/X_old")
    download_command = ["aws", "s3", "sync", path, f"{local_path}/X_old"]
    # Let errors fail the pipeline
    subprocess.run(download_command, check=True)

    logger.info(f"Downloading {obs_path} to {local_path}/old_obs")
    download_command = ["aws", "s3", "sync", obs_path, f"{local_path}/old_obs"]
    # Let errors fail the pipeline
    subprocess.run(download_command, check=True)

    evolve_obs(local_path)
    upload_command = ["aws", "s3", "sync", "--delete", f"{local_path}/new_obs", obs_path]
    subprocess.run(upload_command, check=True)

    params = {
        "kind": "auto",
        "obs_extent": 256,
        "var_extent": 256,
        "cell_order": "row",
        "tile_order": "row",
        "capacity": 1024000,
        "compression": 22,
        "target_array": "X_new",
        "source_array": "X_old",
    }

    executed = evolve_X(cxg=local_path, **params)  # executed is true if a sparse array was upgraded
    if executed:
        logger.info(f"Dataset at {path} computed successfully")
    else:
        logger.info("Dataset was dense, nothing to be done")

    if not dry_run and executed:
        for suffix in ["r", "c"]:
            upload_command = ["aws", "s3", "sync", "--delete", f"{local_path}/X_new" + suffix, path + suffix]
            subprocess.run(upload_command, check=True)

    # Cleanup
    logger.info("Cleaning up local files")
    import shutil

    shutil.rmtree(f"{local_path}/X_old")
    if executed:
        shutil.rmtree(f"{local_path}/X_new")

    shutil.rmtree(f"{local_path}/old_obs")
    shutil.rmtree(f"{local_path}/new_obs")


def evolve_obs(cxg, array_name="old_obs"):
    """
    Creates a new obs array with the new schema
    :param X: the old obs tiledb array
    """
    with tiledb.open(f"{cxg}/{array_name}", "r") as X:
        # load cxg schema
        schema = json.loads(X.meta["cxg_schema"])

        # get all data in the old obs array
        data = X.query().multi_index[:]

        # adjust the schema to expect codes where applicable and store code-to-value mapping
        # dictionary in type_hint["categories"]
        tdb_attrs = []
        new_data = {}
        for attr in X.schema:
            a = attr.name
            type_hint = schema.get(a, {})
            if "categories" in type_hint and len(type_hint.get("categories", [])) > 0.75 * X.shape[0]:
                schema[a]["type"] = "string"
                del schema[a]["categories"]
                tdb_attrs.append(attr)
                new_data[a] = data[a]
            elif "categories" in type_hint:
                cat = pd.Categorical(data[a])
                codes = cat.codes
                new_data[a] = codes
                categories = cat.categories
                schema[a]["categories"] = list(categories)

                dtype = str(cat.codes.dtype)
                tdb_attrs.append(tiledb.Attr(name=a, dtype=dtype, filters=attr.filters))
            else:
                tdb_attrs.append(attr)
                new_data[a] = data[a]

        new_schema = tiledb.ArraySchema(domain=X.schema.domain, attrs=tdb_attrs)
        tiledb.Array.create(
            f"{cxg}/new_obs",
            new_schema,
        )
        with tiledb.open(f"{cxg}/new_obs", "w") as new_X:
            new_X[:] = new_data
            new_X.meta["cxg_schema"] = json.dumps(schema)


def evolve_X(**kwargs):
    """
    Computes the evolved cxg from the sparse source_array and saves it to the target_array destination
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

    st = time.time()
    with tiledb.scope_ctx(ctx):
        with tiledb.open(f"{cxg}/{source_array}", "r") as old_X:
            in_sparse = old_X.schema.sparse
            if not in_sparse:  # exit early, nothing to be done
                return False

            create_new_X(
                cxg,
                target_array,
                compression,
                X_extent,
                old_X.schema,
                cell_order,
                tile_order,
                capacity,
            )
            logger.info("created, starting to read...")
            with tiledb.open(f"{cxg}/{target_array}", "w") as new_X_r:
                with tiledb.open(f"{cxg}/{target_array}c", "w") as new_X_c:
                    i, chunk = 0, 200_000
                    while i < old_X.shape[0]:
                        logger.info(f"Chunk {i}")
                        dat = old_X[i : i + chunk]
                        logger.info("Got dat")
                        new_X_r[dat["obs"]] = {"var": dat["var"], "": dat[""]}
                        new_X_c[dat["var"]] = {"obs": dat["obs"], "": dat[""]}
                        i += chunk

            logger.info("consolidating...")
            for suffix in ["r", "c"]:
                tiledb.consolidate(f"{cxg}/{target_array}{suffix}")
                tiledb.vacuum(f"{cxg}/{target_array}{suffix}")

        logger.info(f"Took {time.time() - st} s")

    logger.info("done.")
    return True


def create_new_X(cxg, target_array, compression, X_extent, old_schema, cell_order, tile_order, capacity):
    logger.info(f"create_new_X: {True}")
    attr_filters = tiledb.FilterList([tiledb.ZstdFilter(level=compression)])
    dim_filters = tiledb.FilterList([tiledb.ByteShuffleFilter(), tiledb.ZstdFilter(level=compression)])
    old_dims = [old_schema.domain.dim(d) for d in range(old_schema.domain.ndim)]
    old_attr = old_schema.attr(0)

    for it, suffix, name_dim, name_attr in zip([0, 1], ["r", "c"], ["obs", "var"], ["var", "obs"]):
        tiledb.Array.create(
            f"{cxg}/{target_array}{suffix}",
            tiledb.ArraySchema(
                domain=tiledb.Domain(
                    [
                        tiledb.Dim(
                            name=name_dim,
                            domain=old_dims[it].domain,
                            tile=X_extent[it],
                            dtype=old_dims[it].dtype,
                            filters=dim_filters,
                        )
                    ]
                ),
                sparse=True,
                allows_duplicates=True,
                attrs=[
                    tiledb.Attr(name="", dtype=old_attr.dtype, filters=attr_filters),
                    tiledb.Attr(name=name_attr, dtype=old_dims[1 - it].dtype, filters=dim_filters),
                ],
                cell_order=f"{cell_order}-major",
                tile_order=f"{tile_order}-major",
                capacity=capacity,
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
