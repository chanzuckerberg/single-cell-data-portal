from typing import List
from functools import cache

import tiledb


@cache
def get_all_dataset_ids(tdb_group: str) -> List[str]:
    with tiledb.open(f"{tdb_group}/obs", "r") as obs:
        all_dataset_ids = obs.query(attrs=[], dims=["dataset_id"]).df[:].dataset_id.unique()
    all_dataset_ids.sort()
    return all_dataset_ids

