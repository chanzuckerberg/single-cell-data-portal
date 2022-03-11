import os

import tiledb

from backend.corpora.common.utils.secret_config import SecretConfig
from backend.wmg.data.utils import frac_mem


class WmgConfig(SecretConfig):
    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="wmg_config", **kwargs)

    def get_defaults_template(self):
        deployment_stage = os.getenv("DEPLOYMENT_STAGE", "test")
        defaults_template = {"bucket": f"wmg-{deployment_stage}"}
        return defaults_template


def create_fast_ctx(config_overrides: dict = {}) -> tiledb.Ctx:
    return create_ctx(fast_config(config_overrides))


def create_ctx(config: dict) -> tiledb.Ctx:
    if not config:
        config = base_config()
    cfg = tiledb.Config(config)
    ctx = tiledb.Ctx(config=cfg)
    return ctx


def base_config():
    return {
        "vfs.s3.region": "us-west-2",
    }


def fast_config(config_overrides: dict = {}) -> dict:
    # consolidation buffer heuristic to prevent thrashing: total_mem/io_concurrency_level, rounded to GB
    io_concurrency_level = int(tiledb.Config()["sm.io_concurrency_level"])
    consolidation_buffer_size = (
            (int(frac_mem(0.1) / io_concurrency_level) + (1024 ** 3 - 1)) // (1024 ** 3) * (1024 ** 3)
    )

    config = {
        "py.init_buffer_bytes": 16 * 1024 ** 3,  # needs to be at least 8GB
        "sm.tile_cache_size": frac_mem(0.5),
        "sm.consolidation.buffer_size": consolidation_buffer_size,
        "sm.query.sparse_unordered_with_dups.non_overlapping_ranges": "true",
        "vfs.s3.region": "us-west-2"
    }
    config.update(config_overrides)
    return config
