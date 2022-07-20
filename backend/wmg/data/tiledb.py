import os

import psutil
import tiledb

from backend.common.utils.math_utils import MB, GB


def create_ctx(config_overrides: dict = {}) -> tiledb.Ctx:
    cfg = {
        "py.init_buffer_bytes": int(0.5 * GB),
        "sm.tile_cache_size": 100 * MB if os.getenv("DEPLOYMENT_STAGE", "test") == "test" else virtual_memory_size(0.5),
        "sm.consolidation.buffer_size": consolidation_buffer_size(0.1),
        "sm.query.sparse_unordered_with_dups.non_overlapping_ranges": "true",
    }
    if boto_endpoint_url := os.getenv("BOTO_ENDPOINT_URL"):
        cfg.update(
            {
                "vfs.s3.endpoint_override": boto_endpoint_url,
                # localstack does not support S3 virtual addressing (per
                # https://docs.aws.amazon.com/AmazonS3/latest/userguide/VirtualHosting.html)
                "vfs.s3.use_virtual_addressing": "false",
            }
        )
    else:
        cfg.update(
            {
                "vfs.s3.region": "us-west-2",
            }
        )
    cfg.update(config_overrides)
    return tiledb.Ctx(config=tiledb.Config(cfg))


def virtual_memory_size(vm_fraction: float) -> int:
    mem_size_bytes = psutil.virtual_memory().total
    fractional_mem_bytes = int(vm_fraction * mem_size_bytes)
    return fractional_mem_bytes // MB * MB  # round down to MB boundary


def consolidation_buffer_size(vm_fraction: float) -> int:
    # consolidation buffer heuristic to prevent thrashing: total_mem/io_concurrency_level, rounded to GB
    io_concurrency_level = int(tiledb.Config()["sm.io_concurrency_level"])
    buffer_size = int(virtual_memory_size(vm_fraction) / io_concurrency_level) + (GB - 1)
    return buffer_size // GB * GB  # round down to GB boundary
