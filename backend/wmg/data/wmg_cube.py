from backend.corpora.common.utils.math_utils import MB
from backend.wmg.data.schemas.corpus_schema import INTEGRATED_ARRAY_NAME, OBS_ARRAY_NAME, VAR_ARRAY_NAME
from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME, EXPRESSION_SUMMARY_CUBE_NAME

import concurrent
import logging

import pandas as pd
import tiledb

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def make_cube_index(tdb_group, cube_dims):
    """
    Create index for queryable dimensions
    """
    with tiledb.open(f"{tdb_group}/{OBS_ARRAY_NAME}") as obs:
        cell_labels = obs.query(use_arrow=False).df[:]
    cell_labels.sort_values(by=["obs_idx"], inplace=True, ignore_index=True)

    cell_labels = pd.DataFrame(
        data={k: cell_labels[k].astype("category") for k in cube_dims},
        index=cell_labels.obs_idx,
    )

    cube_index = pd.DataFrame(cell_labels.value_counts(), columns=["n"])
    cube_index["cube_idx"] = range(len(cube_index))

    cell_labels = cell_labels.join(cube_index.cube_idx, on=cube_dims)

    # we failed to correctly create the cube if these are false
    assert len(cell_labels.index) == cell_labels.index[-1] + 1
    assert cell_labels.index[0] == 0

    return cell_labels, cube_index



