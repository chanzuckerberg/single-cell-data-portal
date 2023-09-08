import os
import sys

from click import Context

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import DatasetId


def backfill_primary_cell_count(ctx: Context, mapping: dict[str, int]) -> None:
    business_logic: BusinessLogic = ctx.obj["business_logic"]
    for i, key in enumerate(mapping):
        dv = business_logic.get_dataset_version_from_canonical(DatasetId(key))
        dv.metadata.primary_cell_count = mapping.get(str(dv.dataset_id), dv.metadata.primary_cell_count)
        print(key, dv.metadata.primary_cell_count, f"({i+1}/{len(mapping)})")
        business_logic.set_dataset_metadata(dv.dataset_id, dv.metadata)
