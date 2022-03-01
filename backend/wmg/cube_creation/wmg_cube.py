"""
TODO:
1. parallelize cube building?  See https://more-itertools.readthedocs.io/en/stable/api.html#more_itertools.chunked
2. calculate bins & mean more efficiently in cube builder
3. large cube tests

"""
from contextlib import contextmanager
import tempfile
from typing import Tuple, Union, Dict
import time
from dataclasses import dataclass
import concurrent.futures
import gc

import tiledb
import numpy as np
import pandas as pd
import fsspec

from .utils import create_fast_ctx, logmsg, create_ctx
from .compute import coo_cube_pass1_into, coo_cube_pass2_into
from .benchmark_consts import feature_ids_most_dense, feature_ids_least_dense, parent_cell_types


big_cube_attrs = [
    "dataset_id",
    "assay_ontology_term_id",
    "development_stage_ontology_term_id",
    "disease_ontology_term_id",
    "ethnicity_ontology_term_id",
    "sex_ontology_term_id",
]
big_cube_dims = [
    "cell_type_ontology_term_id",
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
    *big_cube_attrs,
]
small_cube_dims = [
    "cell_type_ontology_term_id",
    "tissue_ontology_term_id",
    "organism_ontology_term_id",
]


@dataclass
class CubeQuery:
    """
    Very tightly coupled to the cube schema.
    """

    # Dimensions
    feature_id: Union[slice, bytes, list] = slice(None)
    cell_type_ontology_term_id: Union[slice, bytes, list] = slice(None)
    tissue_ontology_term_id: Union[slice, bytes, list] = slice(None)
    organism_ontology_term_id: Union[slice, bytes, list] = slice(None)

    # Attributes (big cube only)
    filters: list[str] = None

    def query(self):
        """return tuple that can be used as a cube query."""
        return (
            self.feature_id,
            self.cell_type_ontology_term_id,
            self.tissue_ontology_term_id,
            self.organism_ontology_term_id,
        )

    def attr_cond(self):
        if self.filters:
            return tiledb.QueryCondition(" and ".join(self.filters))
        else:
            return None


def create_empty_cube(uri: str, other_attrs=[], verbose: bool = False):
    if tiledb.array_exists(uri):
        tiledb.remove(uri)

    filters = [tiledb.ZstdFilter(level=-22)]
    domain = tiledb.Domain(
        [
            tiledb.Dim(name="feature_id", domain=None, tile=None, dtype="ascii", filters=filters),
            tiledb.Dim(name="cell_type_ontology_term_id", domain=None, tile=None, dtype="ascii", filters=filters),
            tiledb.Dim(name="tissue_ontology_term_id", domain=None, tile=None, dtype="ascii", filters=filters),
            tiledb.Dim(name="organism_ontology_term_id", domain=None, tile=None, dtype="ascii", filters=filters),
        ]
    )
    attrs = [
        tiledb.Attr(name="n_cells", dtype=np.uint32, filters=filters),
        tiledb.Attr(name="nnz", dtype=np.uint64, filters=filters),
        tiledb.Attr(name="sum", dtype=np.float32, filters=filters),
        *[tiledb.Attr(name=a, dtype="ascii", var=True, filters=filters) for a in other_attrs],
    ]

    tiledb.Array.create(
        uri,
        tiledb.ArraySchema(
            domain=domain,
            sparse=True,
            allows_duplicates=True,
            attrs=attrs,
            cell_order="row-major",
            tile_order="row-major",
            capacity=10000,
        ),
    )


def make_cube_index(args, cube_dims):
    with tiledb.open(f"{args.tdb_group}/obs") as obs:
        cell_labels = obs.query(use_arrow=False).df[:]
    cell_labels.sort_values(by=['obs_idx'], inplace=True, ignore_index=True)
    cell_labels = pd.DataFrame(
        data={k: cell_labels[k].astype("category") for k in cube_dims},
        index=cell_labels.obs_idx,
    )
    cube_index = pd.DataFrame(cell_labels.value_counts(), columns=["n"])
    cube_index["cube_idx"] = range(len(cube_index))

    cell_labels = cell_labels.join(cube_index.cube_idx, on=cube_dims)

    # we failed to correctly create the corpus if these are false
    assert len(cell_labels.index) == cell_labels.index[-1] + 1
    assert cell_labels.index[0] == 0

    return cell_labels, cube_index


def reduce_X(args, start_time, cube_indices, reducer, *accum):
    verbose = args.verbose

    with concurrent.futures.ThreadPoolExecutor() as tp:
        cfg = {
            "py.init_buffer_bytes": 512 * 1024 ** 2,
            "py.exact_init_buffer_bytes": "true",
        }
        with tiledb.open(f"{args.tdb_group}/raw", ctx=create_fast_ctx(cfg)) as X:
            iterable = X.query(return_incomplete=True, order="U", attrs=["data"])
            future = None
            for i, result in enumerate(iterable.df[:]):
                if verbose:
                    logmsg(f"reduce raw X, iter {i}, {time.time()-start_time}")
                if future is not None:
                    future.result()  # forces a wait
                future = tp.submit(
                    reducer,
                    result["data"].values,
                    result["obs_idx"].values,
                    result["var_idx"].values,
                    cube_indices,
                    *accum,
                )

        return accum


def build_in_mem_cube(feature_ids, cube_index, other_attrs, cube_sum, cube_nnz, verbose=False):
    if verbose:
        logmsg("Building in-mem cube")

    # Count total values so we can allocate buffers once
    total_vals = 0
    for cube_idx in cube_index.cube_idx.values:
        mask = cube_nnz[cube_idx] != 0
        total_vals += np.count_nonzero(mask)

    # allocate buffers
    dims = [
        np.empty((total_vals,), dtype=object),
        np.empty((total_vals,), dtype=object),
        np.empty((total_vals,), dtype=object),
        np.empty((total_vals,), dtype=object),
    ]
    vals = {
        "sum": np.empty((total_vals,)),
        "nnz": np.empty((total_vals,), dtype=np.uint64),
        "n_cells": np.empty((total_vals,), dtype=np.uint32),
        **{k: np.empty((total_vals,), dtype=object) for k in other_attrs},
    }

    # populate buffers
    idx = 0
    for grp in cube_index.to_records():
        (
            cell_type_ontology_term_id,
            tissue_ontology_term_id,
            organism_ontology_term_id,
            *attr_values,
            n,
            cube_idx,
        ) = grp.tolist()
        mask = cube_nnz[cube_idx] != 0
        n_vals = np.count_nonzero(mask)
        if n_vals == 0:
            continue

        if verbose > 2:
            logmsg(grp)

        dims[0][idx : idx + n_vals] = feature_ids.feature_id.values[mask]
        dims[1][idx : idx + n_vals] = cell_type_ontology_term_id
        dims[2][idx : idx + n_vals] = tissue_ontology_term_id
        dims[3][idx : idx + n_vals] = organism_ontology_term_id

        vals["sum"][idx : idx + n_vals] = cube_sum[cube_idx, mask]
        vals["nnz"][idx : idx + n_vals] = cube_nnz[cube_idx, mask]
        vals["n_cells"][idx : idx + n_vals] = n  # wasteful

        for i, k in enumerate(other_attrs):
            vals[k][idx : idx + n_vals] = attr_values[i]

        idx += n_vals

    return dims, vals


def load_small_cube(args, uri: str):
    verbose = args.verbose
    ctx = create_fast_ctx()
    start_time = time.time()

    if verbose:
        logmsg("start")

    with tiledb.open(f"{args.tdb_group}/var", ctx=ctx) as var:
        feature_ids = var.query(dims=["feature_id"], attrs=[]).df[:]
    n_features = len(feature_ids)

    ##
    ## Reduce X
    ##
    cell_labels, cube_index = make_cube_index(args, small_cube_dims)
    n_groups = len(cube_index)

    cube_sum = np.zeros((n_groups, n_features), dtype=np.float32)
    cube_nnz = np.zeros((n_groups, n_features), dtype=np.uint64)
    cube_min = np.zeros((n_groups, n_features), dtype=np.float32)
    cube_max = np.zeros((n_groups, n_features), dtype=np.float32)

    # pass 1 - sum, nnz, min, max
    reduce_X(args, start_time, cell_labels.cube_idx.values, coo_cube_pass1_into, cube_sum, cube_nnz, cube_min, cube_max)

    return build_in_mem_cube(feature_ids, cube_index, [], cube_sum, cube_nnz, verbose)


def load_big_cube(args, uri: str):
    verbose = args.verbose
    ctx = create_fast_ctx()
    start_time = time.time()

    if verbose:
        logmsg("start")

    with tiledb.open(f"{args.tdb_group}/var", ctx=ctx) as var:
        feature_ids = var.query(dims=["feature_id"], attrs=[], use_arrow=False).df[:]
    n_features = len(feature_ids)

    ##
    ## Reduce X
    ##
    cell_labels, cube_index = make_cube_index(args, big_cube_dims)
    n_groups = len(cube_index)

    cube_sum = np.zeros((n_groups, n_features), dtype=np.float32)
    cube_nnz = np.zeros((n_groups, n_features), dtype=np.uint64)
    cube_min = np.zeros((n_groups, n_features), dtype=np.float32)
    cube_max = np.zeros((n_groups, n_features), dtype=np.float32)

    # pass 1 - sum, nnz, min, max
    reduce_X(args, start_time, cell_labels.cube_idx.values, coo_cube_pass1_into, cube_sum, cube_nnz, cube_min, cube_max)

    return build_in_mem_cube(feature_ids, cube_index, big_cube_attrs, cube_sum, cube_nnz, verbose)


def create_small_cube(args):
    verbose = args.verbose
    uri = f"{args.tdb_group}/small-cube"
    start_time = time.time()

    ctx = create_fast_ctx()
    with tiledb.scope_ctx(ctx):

        # Create cube
        create_empty_cube(uri, [], verbose=verbose)

        # load data
        dims, vals = load_small_cube(args, uri)

        if verbose:
            logmsg("Saving small cube to tiledb")
        with tiledb.open(uri, "w") as cube:
            cube[tuple(dims)] = vals

        if verbose:
            logmsg("small cube created, start consolidation")
        tiledb.consolidate(uri)

        if verbose:
            logmsg("small cube consolidated, start vacuumming")
        tiledb.vacuum(uri)

    if verbose:
        logmsg("small cube complete")
    create_cube_sec = time.time() - start_time
    logmsg(f"Small cube: time to create {create_cube_sec}, uri={uri}")


def create_big_cube(args):
    verbose = args.verbose
    uri = f"{args.tdb_group}/big-cube"
    start_time = time.time()

    with tiledb.scope_ctx(create_fast_ctx()):

        # Create cube
        create_empty_cube(uri, big_cube_attrs, verbose=verbose)

        # load data
        dims, vals = load_big_cube(args, uri)

        if verbose:
            logmsg("Saving big cube to tiledb")
        with tiledb.open(uri, "w") as cube:
            try:
                logmsg("why wont this work?")
                cube[tuple(dims)] = vals
            except Exception as e:
                logmsg(e)

        if verbose:
            logmsg("big cube created, start consolidation")
        tiledb.consolidate(uri)

        if verbose:
            logmsg("big cube consolidated, start vacuumming")
        tiledb.vacuum(uri)

    if verbose:
        logmsg("big cube complete")
    create_cube_sec = time.time() - start_time
    logmsg(f"Big cube: time to create {create_cube_sec}, uri={uri}")


def create_test_matrix(rng, args, is_big_cube=False):
    """
    spec: https://docs.google.com/document/d/1xA5ubOqVqZVDlOjtptVkvxWltiY-UogQbBLIYYLAfrE/edit#heading=h.8yyau2t8pqdb

    Filter to a set of cells based on the following cell-label filters:
        * single_tissue_a: tissue_ontology_term_id == "UBERON:0000947"
        * single_tissue_b: tissue_ontology_term_id == "UBERON:0000178"
        * multi_tissue: tissue_ontology_term_id in ["UBERON:0016540", "UBERON:0002894", "UBERON:0004167", "UBERON:0016530", "UBERON:0004725", "UBERON:0000451"]
        * organism_and_tissue: organism_ontology_term_id == "NCBITaxon:10090" AND tissue_ontology_term_id == "UBERON:0000947"
        * all (i.e. no filter)
    cell_type_ontology_term_id cases:
        * individual: Every individual cell type present in the data (ie, unique labels on cells), queried individually but reported as min, max, deciles
        * parent groups: Parent cell types to represent hierarchical-like results:
            * misc_neurons: [CL:0000679, CL:0000748, CL:0000540]
            * misc_macrophage: [CL:0000235, CL:1001603, CL:0002150, CL:0000864, CL:1001603].
            * misc_epithelial: [CL:0000066, CL:0000067, CL:0000068, CL:0000079 ].
        * all: All cell types in the data (i.e. no filter on cell type)
    Filter to a set of genes based on:
        * random: Randomly 50, 100, 200, 500 (separate trials, repeat 5x and report min, max, median)
        * densest: Specifically, the top (1, 10, 30, 100) densest genes (see spreadsheet below)
        * sparsest: Specifically, the top (1, 10, 30, 100) sparsest genes (see spreadsheet below)

    Calculate, separately:
        * Sum of expression (X)
        * NNZ
        * N_cells

    In summary, each test is:
        * Input: 1. a filter, 2. cell cell_type_ontology_term_id(s), 3. gene_name(s) .
        * Output: matrix of cell_type_ontology_term_id, gene_name, value (mean or histogram)
    Also include the count of non-zero non-nan expression values (post-filtering, pre-reduction, i.e. the same count that you would need to compute the mean)

    """
    tdb_group = args.tdb_group
    quick = args.quick

    with tiledb.open(f"{tdb_group}/var") as var:
        features = var.query(dims=["feature_id"], attrs=[]).df[:]
        all_feature_ids = features.feature_id.unique()
    all_feature_ids = np.char.encode(all_feature_ids.astype(str), encoding="ascii")

    with tiledb.open(f"{tdb_group}/obs") as obs:
        labels = obs.query(dims=["cell_type_ontology_term_id"], attrs=[]).df[:]
        all_cell_types = labels.cell_type_ontology_term_id.unique().astype(str)
        all_cell_types = np.char.encode(all_cell_types, encoding="ascii")

    label_filters = {
        "all": {},
        "single_tissue_a": dict(tissue_ontology_term_id=b"UBERON:0000947"),
        "single_tissue_b": dict(tissue_ontology_term_id=b"UBERON:0000178"),
        "multi_tissue": dict(
            tissue_ontology_term_id=[
                b"UBERON:0016540",
                b"UBERON:0002894",
                b"UBERON:0004167",
                b"UBERON:0016530",
                b"UBERON:0004725",
                b"UBERON:0000451",
            ]
        ),
        "organism_and_tissue": dict(
            tissue_ontology_term_id=b"UBERON:0000947",
            organism_ontology_term_id=b"NCBITaxon:10090",
        ),
    }

    def feature_id_filters(trials=1):
        if not quick:
            for n in [1, 10, 30, 100]:
                most_dense = feature_ids_most_dense[-n:]
                yield f"densest_{n}", dict(feature_id=most_dense)

                least_dense = feature_ids_least_dense[0:n]
                yield f"sparsest_{n}", dict(feature_id=least_dense)

        for n in [200, 500] if quick else [50, 100, 200, 500]:
            for r in range(trials):
                feature_choices = rng.choice(len(all_feature_ids), n)
                ids = all_feature_ids[feature_choices].tolist()
                yield f"random_{n}", dict(feature_id=ids)

    def matrix():
        """
        returns:
            filter name	- name of cell filter, str
            cell type name - name of cell type test, str
            genes name - name of feature filter, str
            trials query generator
        """
        # # grouped/parent cell types
        for label_filter_name, label_filter in label_filters.items():
            for cell_type_name, cell_types in parent_cell_types.items():
                for genes_filter_name, genes_filter in feature_id_filters(trials=5):
                    yield label_filter_name, cell_type_name, genes_filter_name, CubeQuery(
                        **(label_filter | genes_filter | dict(cell_type_ontology_term_id=cell_types))
                    )

        # individual cell types aggregated
        if not quick:
            for label_filter_name, label_filter in label_filters.items():
                for cell_type in all_cell_types:
                    cell_types = dict(cell_type_ontology_term_id=cell_type)
                    for genes_filter_name, genes_filter in feature_id_filters(trials=1):
                        yield label_filter_name, "individual", genes_filter_name, CubeQuery(
                            **(label_filter | genes_filter | cell_types)
                        )

    return matrix


def cube_query(cube, query: CubeQuery, verbose: bool = False):
    if verbose:
        print(".", end="", flush=True)
        if verbose >= 3:
            print(query)

    query_index = query.query()
    query_attr_cond = query.attr_cond()

    start_time = time.time()
    qres = cube.query(
        dims=["feature_id", "cell_type_ontology_term_id"],
        attrs=["sum", "n_cells", "nnz"],
        attr_cond=query_attr_cond,
        use_arrow=False,
    ).df[query_index]
    mean_result = qres.groupby(by=["feature_id", "cell_type_ontology_term_id"]).agg(
        {"sum": ("sum"), "nnz": ("sum"), "n_cells": ("sum")}
    )
    mean_t = time.time() - start_time

    return mean_t, mean_result


def push_result(results, filters_name, cell_types_name, gene_name, mean_t, bin_t, nnz, res_count):
    if results is None:
        # initialize
        results = {"filters": [], "cell_types": [], "genes": [], "res_count": [], "nnz": [], "mean_t": []}
    results["filters"].append(filters_name)
    results["cell_types"].append(cell_types_name)
    results["genes"].append(gene_name)
    results["nnz"].append(nnz)
    results["res_count"].append(res_count)
    results["mean_t"].append(mean_t)
    return results


@contextmanager
def preload_cube(args, uri, benchmark_ctx):
    if args.use_mem:
        # copy to /dev/shm
        with tempfile.TemporaryDirectory(dir="/dev/shm/") as tmpd:
            fs = fsspec.open(uri).fs
            if args.verbose >= 2:
                print("preload cube", uri, "to local cache", tmpd)
                start_time = time.time()
            fs.get(uri, f"{tmpd}/", recursive=True)
            if args.verbose >= 2:
                print(f"preload done, took {time.time()-start_time}s")
            yield tmpd

    else:
        if args.verbose >= 2:
            print("preload cube", uri)
            start_time = time.time()
        with tiledb.scope_ctx(create_ctx()):
            with tiledb.open(f"{args.tdb_group}/var") as var:
                features = var.query(dims=["feature_id"], attrs=[]).df[:]
                all_feature_ids = features.feature_id.unique()
            all_feature_ids = np.char.encode(all_feature_ids.astype(str), encoding="ascii")
        with tiledb.scope_ctx(benchmark_ctx):
            with tiledb.open(uri) as cube:
                cube.multi_index[[all_feature_ids[i] for i in range(0, len(all_feature_ids), 40)], :]
        if args.verbose >= 2:
            print(f"preload done, took {time.time()-start_time}s")

        yield uri


def benchmark_cube(args, uri, test_matrix):
    results = None

    ctx = create_ctx({"sm.tile_cache_size": 16 * 1024 ** 3})
    with preload_cube(args, uri, ctx) as cube_uri:
        with tiledb.scope_ctx(ctx):
            with tiledb.open(cube_uri) as cube:
                for filters_name, cell_types_name, gene_name, query in test_matrix():
                    mean_t, mean_result = cube_query(cube, query, args.verbose)
                    nnz = mean_result.nnz.sum()
                    res_count = len(mean_result)
                    results = push_result(results, filters_name, cell_types_name, gene_name, mean_t, 0, nnz, res_count)
                if args.verbose >= 2:
                    print(tiledb.default_ctx().get_stats())

    if args.verbose:
        print()
    df = pd.DataFrame(data=results)
    df = df.groupby(by=["filters", "cell_types", "genes"]).agg(
        res_count_mean=("res_count", "mean"),
        trials=("nnz", "count"),
        nnz_min=("nnz", "min"),
        nnz_max=("nnz", "max"),
        nnz_mean=("nnz", "mean"),
        mean_t_mean=("mean_t", "mean"),
        mean_t_min=("mean_t", "min"),
        mean_t_max=("mean_t", "max"),
        mean_t_med=("mean_t", "median"),
        mean_t_p75=("mean_t", lambda x: x.quantile(0.75)),
        mean_t_p90=("mean_t", lambda x: x.quantile(0.9)),
    )
    df.reset_index(inplace=True)
    print(df.to_csv(index=False))


def benchmark_small_cube(args, rng):
    test_matrix = create_test_matrix(rng, args)
    benchmark_cube(args, f"{args.tdb_group}/small-cube", test_matrix)


def benchmark_big_cube(args, rng):
    test_matrix = create_test_matrix(rng, args, True)
    benchmark_cube(args, f"{args.tdb_group}/big-cube", test_matrix)
