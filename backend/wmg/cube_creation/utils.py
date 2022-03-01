from typing import List
from functools import cache
import datetime

import tiledb
import psutil


@cache
def get_all_dataset_ids(tdb_group: str) -> List[str]:
    with tiledb.open(f"{tdb_group}/obs", "r") as obs:
        all_dataset_ids = obs.query(attrs=[], dims=["dataset_id"]).df[:].dataset_id.unique()
    all_dataset_ids.sort()
    return all_dataset_ids


def base_config():
    return {
        "vfs.s3.region": "us-west-2",
    }


def create_ctx(config: dict = {}) -> tiledb.Ctx:
    cfg = tiledb.Config(base_config() | config)
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

    config = base_config() | {
        "py.init_buffer_bytes": 16 * 1024 ** 3,  # needs to be at least 8GB
        "sm.tile_cache_size": frac_mem(0.5),
        "sm.consolidation.buffer_size": consolidation_buffer_size,
        "sm.query.sparse_unordered_with_dups.non_overlapping_ranges": "true",
    }
    config.update(config_overrides)
    return config


def create_fast_ctx(config_overrides: dict = {}) -> tiledb.Ctx:
    return create_ctx(fast_config(config_overrides))


# def logmsg(*args) -> None:
#     fmt = "%Y-%m-%d:%H:%M:%S"
#     print(f"[{datetime.datetime.now().strftime(fmt)}]:", *args)
