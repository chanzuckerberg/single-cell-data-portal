import os
import sys

from click import Context

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.layers.business.business import BusinessLogic
from backend.layers.common.entities import DatasetMetadata, DatasetVersionId


def backfill_primary_cell_count(ctx: Context, mapping: dict[str, int]) -> None:
    business_logic: BusinessLogic = ctx.obj["business_logic"]
    dataset_versions = business_logic.get_dataset_versions_by_id([DatasetVersionId(key) for key in map])
    for dv in dataset_versions:
        dataset_metadata = DatasetMetadata(**dv.dataset_metadata)
        dataset_metadata.primary_cell_count = mapping.get(dv.id, dataset_metadata.primary_cell_count)
        business_logic.set_dataset_metadata(dv.id, dataset_metadata)
