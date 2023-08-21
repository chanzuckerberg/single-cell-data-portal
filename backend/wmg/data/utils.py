import logging
import os
import time
from typing import Dict, List

import numpy as np
import requests
import tiledb
import yaml
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

from backend.wmg.data.constants import CL_PINNED_CONFIG_URL, WMG_PINNED_SCHEMA_VERSION
from backend.wmg.data.schemas.corpus_schema import OBS_ARRAY_NAME


def log_func_runtime(func):
    # This decorator function logs the execution time of the function object passed
    def wrap_func(*args, **kwargs):
        logger = logging.getLogger(func.__module__)
        start = time.perf_counter()
        result = func(*args, **kwargs)
        stop = time.perf_counter()
        logger.info(f"Function {func.__name__} executed in {(stop-start):.4f}s")
        return result

    return wrap_func


def get_all_dataset_ids(tdb_group: str) -> List[str]:
    with tiledb.open(f"{tdb_group}/{OBS_ARRAY_NAME}", "r") as obs:
        all_dataset_ids = obs.query(attrs=[], dims=["dataset_id"]).df[:].dataset_id.unique()
    all_dataset_ids.sort()
    return all_dataset_ids


@log_func_runtime
def get_expression_summary_cube_gene_count(tbd_group: str) -> int:
    with tiledb.open(tbd_group) as obs:
        gene_count = len(obs.query(dims=["gene_ontology_term_id"]).df[:].gene_ontology_term_id.unique())
    return gene_count


@log_func_runtime
def get_cell_count_cube_count(tbd_group: str) -> int:
    with tiledb.open(tbd_group) as obs:
        cell_count = obs.query(attrs=["n_cells"]).df[:].n_cells.sum()
    return cell_count


def create_empty_cube(uri: str, schema):
    """
    Create an empty cube with expected schema (dimensions and attributes) at given uri
    """
    tiledb.Array.create(uri, schema, overwrite=True)


def find_all_dim_option_values(snapshot, organism: str, dimension: str) -> list:
    all_filter_options = set()
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


def get_datasets_from_curation_api():
    # hardcode to staging backend if deployment is rdev or test
    API_URL = (
        "https://api.cellxgene.staging.single-cell.czi.technology"
        if os.environ.get("DEPLOYMENT_STAGE") in ["test", "rdev"]
        else os.getenv("API_URL")
    )

    datasets = {}
    if API_URL:
        session = setup_retry_session()
        dataset_metadata_url = f"{API_URL}/curation/v1/datasets?schema_version={WMG_PINNED_SCHEMA_VERSION}"
        response = session.get(dataset_metadata_url)
        if response.status_code == 200:
            datasets = response.json()
    return datasets


def get_collections_from_curation_api():
    # hardcode to staging backend if deployment is rdev or test
    API_URL = (
        "https://api.cellxgene.staging.single-cell.czi.technology"
        if os.environ.get("DEPLOYMENT_STAGE") in ["test", "rdev"]
        else os.getenv("API_URL")
    )

    collections = {}
    if API_URL:
        session = setup_retry_session()
        dataset_metadata_url = f"{API_URL}/curation/v1/collections"
        response = session.get(dataset_metadata_url)
        if response.status_code == 200:
            collections = response.json()
    return collections


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
    d = dict(zip(np.unique(a), [list(set(x)) for x in slists]))
    return d


def get_pinned_ontology_url(name: str):
    """
    This function retrieves the URL of the pinned ontology based on the provided name.

    Parameters:
    name (str): The name of the ontology (e.g. cl-basic.obo).

    Returns:
    str: The URL of the pinned ontology.
    """
    session = setup_retry_session()
    response = session.get(CL_PINNED_CONFIG_URL)
    response.raise_for_status()
    decoded_yaml = yaml.safe_load(response.content.decode())
    key = decoded_yaml["CL"]["latest"]
    cl_url = decoded_yaml["CL"]["urls"][key]
    cl_url = cl_url.split("cl.owl")[0] + name
    return cl_url
