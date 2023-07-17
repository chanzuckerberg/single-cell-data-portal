import os
from typing import Optional

import tiledb

from backend.common.utils.math_utils import GB, MB
from backend.common.utils.tiledb import consolidation_buffer_size, virtual_memory_size


def create_ctx(config_overrides: Optional[dict] = None) -> tiledb.Ctx:
    config_overrides = config_overrides or {}
    cfg = {
        "py.init_buffer_bytes": int(0.05 * GB),
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
