import os
from functools import lru_cache
from typing import Dict

import numpy as np
import pandas as pd
import requests
from cellxgene_ontology_guide.ontology_parser import OntologyParser
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

from backend.common.census_cube.data.constants import (
    CENSUS_CUBE_API_SNAPSHOT_FS_CACHE_ROOT_PATH,
    CENSUS_CUBE_DATA_SCHEMA_VERSION,
    CENSUS_CUBE_PINNED_SCHEMA_VERSION,
)
from backend.common.census_cube.data.snapshot import load_snapshot
from backend.common.constants import DEPLOYMENT_STAGE_TO_API_URL
from backend.common.utils.rollup import rollup_across_cell_type_descendants

DEPLOYMENT_STAGE = os.environ.get("DEPLOYMENT_STAGE", "")
SNAPSHOT_FS_ROOT_PATH = CENSUS_CUBE_API_SNAPSHOT_FS_CACHE_ROOT_PATH if (DEPLOYMENT_STAGE != "test") else None

# exported and used by all modules related to the census cube
ontology_parser = OntologyParser(schema_version=f"v{CENSUS_CUBE_PINNED_SCHEMA_VERSION}")


def find_all_dim_option_values(snapshot, organism: str, dimension: str) -> list:
    organism_key = "organism_ontology_term_id__" + organism
    all_filter_options = snapshot.filter_relationships[organism_key].get(dimension, [])
    return [option.split("__")[1] for option in all_filter_options]


def find_dim_option_values(criteria: Dict, snapshot, dimension: str) -> list:
    """Find values for the specified dimension that satisfy the given filtering criteria,
    ignoring any criteria specified for the given dimension."""

    filter_options_criteria = dict(criteria)
    # Remove gene_ontology_term_ids from the criteria as it is not an eligible cross-filter dimension.
    filter_options_criteria.pop("gene_ontology_term_ids", None)

    # depluralize `dimension` if necessary
    dimension = depluralize(dimension)

    # each element  in `linked_filter_sets` corresponds to the set of filters linked to the attributes specified for a corresponding criteria key
    linked_filter_sets = []

    # `all_criteria_attributes` is the set of all attributes specified across all criteria
    all_criteria_attributes = set()

    for key in filter_options_criteria:
        attrs = filter_options_criteria[key]

        # depluralize `key` if necessary
        key = depluralize(key)

        # ignore the criteria for the specified dimension
        if key != dimension:
            if isinstance(attrs, list):
                if len(attrs) > 0:
                    # prepend the key to each attribute value
                    prefixed_attributes = [key + "__" + val for val in attrs]
                    all_criteria_attributes = all_criteria_attributes.union(prefixed_attributes)

                    # for each attribute (attr) in `prefixed_attributes`,
                    # get the set of filters for the specified dimension that are linked to `attr`
                    linked_filter_set = set()
                    for attr in prefixed_attributes:
                        if dimension in snapshot.filter_relationships.get(attr, {}):
                            linked_filter_set = linked_filter_set.union(
                                set(snapshot.filter_relationships[attr][dimension])
                            )

                    linked_filter_sets.append(linked_filter_set)
            else:
                if attrs != "":
                    prefixed_attribute = key + "__" + attrs
                    all_criteria_attributes.add(prefixed_attribute)
                    if dimension in snapshot.filter_relationships.get(prefixed_attribute, {}):
                        linked_filter_sets.append(set(snapshot.filter_relationships[prefixed_attribute][dimension]))

    # the candidate options are the intersection of the sets of linked filters for each criteria key
    if len(linked_filter_sets) > 1:
        candidate_options = linked_filter_sets[0].intersection(*linked_filter_sets[1:])
    else:
        candidate_options = linked_filter_sets[0]

    # each valid option MUST be linked to at least one attribute specified in the criteria
    # otherwise, there will be no data to display if that particular option is selected because
    # the intersection will be null.
    valid_options = []
    for v in candidate_options:
        loop_back_options = snapshot.filter_relationships.get(v, {})
        all_loop_back_options = []
        for dim in loop_back_options:
            all_loop_back_options.extend(loop_back_options[dim])

        if len(set(all_loop_back_options).intersection(all_criteria_attributes)) > 0:
            valid_options.append(v)

    # remove the prefix from each valid option and return the result
    return [option.split("__")[1] for option in valid_options]


def depluralize(x):
    return x[:-1] if x[-1] == "s" else x


def setup_retry_session(retries=3, backoff_factor=2, status_forcelist=(500, 502, 503, 504), method_whitelist=None):
    session = requests.Session()

    if method_whitelist is None:
        method_whitelist = {"GET"}

    retry = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
        allowed_methods=method_whitelist,
    )

    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)

    return session


def get_datasets_from_discover_api():
    # hardcode to staging backend if deployment is rdev or test
    deployment_stage = os.environ.get("DEPLOYMENT_STAGE")
    API_URL = DEPLOYMENT_STAGE_TO_API_URL.get(
        deployment_stage, "https://api.cellxgene.staging.single-cell.czi.technology"
    )

    datasets = {}
    if API_URL:
        session = setup_retry_session()
        dataset_metadata_url = f"{API_URL}/curation/v1/datasets?schema_version={CENSUS_CUBE_PINNED_SCHEMA_VERSION}"
        response = session.get(dataset_metadata_url)
        if response.status_code == 200:
            datasets = response.json()
    return datasets


def get_collections_from_discover_api():
    # hardcode to staging backend if deployment is rdev or test
    deployment_stage = os.environ.get("DEPLOYMENT_STAGE")
    API_URL = DEPLOYMENT_STAGE_TO_API_URL.get(
        deployment_stage, "https://api.cellxgene.staging.single-cell.czi.technology"
    )

    collections = {}
    if API_URL:
        session = setup_retry_session()
        dataset_metadata_url = f"{API_URL}/curation/v1/collections"
        response = session.get(dataset_metadata_url)
        if response.status_code == 200:
            collections = response.json()
    return collections


def build_filter_relationships(cell_counts_df: pd.DataFrame):
    # get a dataframe of the columns that are not numeric
    df_filters = cell_counts_df.select_dtypes(exclude="number")
    # get a numpy array of the column names with shape (1, n_cols)
    cols = df_filters.columns.values[None, :]

    # tile the column names row-wise to match the shape of the dataframe and concatenate to the values
    # this ensures that filter values will never collide across columns.
    mat = np.tile(cols, (cell_counts_df.shape[0], 1)) + "__" + df_filters.values

    # for each cell, get all pairwise combinations of filters compresent in that cell
    # these are the edges of the filter relationships graph
    Xs = []
    Ys = []
    for i in range(mat.shape[0]):
        Xs.extend(np.repeat(mat[i], mat[i].size))
        Ys.extend(np.tile(mat[i], mat[i].size))

    # get all the unique combinations of filters
    Xs, Ys = np.unique(np.array((Xs, Ys)), axis=1)

    # exclude self-relationships
    filt = Xs != Ys
    Xs, Ys = Xs[filt], Ys[filt]

    # convert the edges to a linked list representation
    filter_relationships_linked_list = to_dict(Xs, Ys)

    # reorganize the linked list representation to a nested linked list representation
    # where the filter columns are separated into distinct dictionaries
    # e.g. instead of {"cell_type_ontology_term_id__beta cell": ["dataset_id__Single cell transcriptome analysis of human pancreas", "assay_ontology_term_id__assay_type", ...], ...},
    # it's now {"cell_type_ontology_term_id__beta cell": {"dataset_id": ["dataset_id__Single cell transcriptome analysis of human pancreas", ...], "assay_ontology_term_id": ["assay_ontology_term_id__assay_type", ...], ...}, ...}.
    # This structure is easier to parse by the `/query` endpoint.
    for k, v in filter_relationships_linked_list.items():
        filter_relationships_linked_list[k] = to_dict([x.split("__")[0] for x in v], v)

    return filter_relationships_linked_list


def to_dict(a, b):
    """
    convert a flat key array (a) and a value array (b) into a dictionary with values grouped by keys
    """
    a = np.array(a)
    b = np.array(b)
    idx = np.argsort(a)
    a = a[idx]
    b = b[idx]
    bounds = np.where(a[:-1] != a[1:])[0] + 1
    bounds = np.append(np.append(0, bounds), a.size)
    bounds_left = bounds[:-1]
    bounds_right = bounds[1:]
    slists = [b[bounds_left[i] : bounds_right[i]] for i in range(bounds_left.size)]
    d = dict(zip(np.unique(a), [list(set(x)) for x in slists], strict=False))
    return d


@lru_cache(maxsize=None)
def get_all_cell_type_ids_in_corpus(root_node="CL:0000000") -> list[str]:
    """
    Retrieve all cell type IDs in the corpus that have at least one cell present, starting from a specified root node in the ontology.

    This function uses the ontology parser to get all descendant cell type IDs from the root node, including the root node itself.
    It then fetches the cell counts dataframe from the snapshot and groups it by cell type ontology term ID to sum up the cell counts.
    Cell types with zero cells are added to ensure all cell types from the ontology are represented in the output, even if they have no cells.
    Descendant cell counts are rolled up into each cell type node and remaining cell types with zero cells after rollup are filtered out.
    The function finally returns a list of cell type IDs that have cells.

    Parameters:
        root_node (str): The ontology term ID of the root node from which to start gathering cell type IDs.

    Returns:
        list[str]: A list of cell type ontology term IDs that have at least one cell in the corpus.
    """
    snapshot = load_snapshot(
        snapshot_schema_version=CENSUS_CUBE_DATA_SCHEMA_VERSION,
        explicit_snapshot_id_to_load=None,
        snapshot_fs_root_path=SNAPSHOT_FS_ROOT_PATH,
    )
    all_cell_type_ids = ontology_parser.get_term_descendants(root_node, include_self=True)
    cell_counts_df = snapshot.cell_counts_df

    cell_counts_df = (
        cell_counts_df.groupby("cell_type_ontology_term_id").sum(numeric_only=True)[["n_cells"]].reset_index()
    )
    # to_attach is a DataFrame that will contain cell type ontology term ids that are not present in the input cell counts dataframe. These will be added with 0 counts.
    to_attach = pd.DataFrame()
    to_attach["cell_type_ontology_term_id"] = [
        i for i in all_cell_type_ids if i not in cell_counts_df["cell_type_ontology_term_id"].values
    ]
    to_attach["n_cells"] = 0
    cell_counts_df = pd.concat([cell_counts_df, to_attach], axis=0)
    cell_counts_df_rollup = rollup_across_cell_type_descendants(cell_counts_df)
    cell_counts_df_rollup = cell_counts_df_rollup.set_index("cell_type_ontology_term_id")["n_cells"]
    all_cell_type_ids_in_corpus = cell_counts_df_rollup.index.values[cell_counts_df_rollup.values > 0]
    return list(all_cell_type_ids_in_corpus)


@lru_cache(maxsize=None)
def get_all_tissue_ids_in_corpus() -> list[str]:
    """
    Retrieve all tissue IDs in the corpus that have at least one cell present.

    This function fetches the cell counts dataframe from the snapshot.
    It then uses a utility function to convert the grouped data into a dictionary mapping tissue IDs to cell type IDs.
    The function finally returns a list of tissue IDs that have cells.

    Returns:
        list[str]: A list of tissue ontology term IDs that have at least one cell in the corpus.
    """
    snapshot = load_snapshot(
        snapshot_schema_version=CENSUS_CUBE_DATA_SCHEMA_VERSION,
        explicit_snapshot_id_to_load=None,
        snapshot_fs_root_path=SNAPSHOT_FS_ROOT_PATH,
    )
    cell_counts_df = snapshot.cell_counts_df
    uberon_by_celltype = to_dict(
        cell_counts_df["tissue_ontology_term_id"], cell_counts_df["cell_type_ontology_term_id"]
    )
    return list(uberon_by_celltype.keys())
