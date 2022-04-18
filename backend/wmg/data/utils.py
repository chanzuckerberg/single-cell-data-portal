from functools import cache
from typing import List

import tiledb

from backend.wmg.data.schemas.corpus_schema import OBS_ARRAY_NAME


@cache
def get_all_dataset_ids(tdb_group: str) -> List[str]:
    with tiledb.open(f"{tdb_group}/{OBS_ARRAY_NAME}", "r") as obs:
        all_dataset_ids = obs.query(attrs=[], dims=["dataset_id"]).df[:].dataset_id.unique()
    all_dataset_ids.sort()
    return all_dataset_ids
