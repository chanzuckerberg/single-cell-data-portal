import argparse
import io
import os
import time
import random
import pprint
from cProfile import Profile
from pstats import Stats

import tiledb
from tiledb import Array
import pandas as pd

from backend.wmg.config import WmgConfig
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.api.v2 import *

def local_disk_snapshot():
    return WmgSnapshot(snapshot_identifier="",
                       expression_summary_cube=tiledb.open('staging-snapshot/expression_summary'),
                       expression_summary_default_cube=tiledb.open('staging-snapshot/expression_summary_default'))

def s3_snapshot():
    def tiledb_open_s3_uri(s3_uri):
        # TODO (fix): Is this config correct?
        return tiledb.open(
            s3_uri,
            config=tiledb.Config(
                {"vfs.s3.region": "us-west-2", 'py.init_buffer_bytes': 128*1024**2}))
    
    return WmgSnapshot(snapshot_identifier="",
                       expression_summary_cube=tiledb_open_s3_uri('s3://cellxgene-wmg-staging/snapshots/v3/1701021778/expression_summary'),
                       expression_summary_default_cube=tiledb_open_s3_uri('s3://cellxgene-wmg-staging/snapshots/v3/1701021778/expression_summary_default'))

def print_profile_results_web_app_query(snapshot, profile_func, randomize_reqs, num_trials, wait):
    expr_default_cube = snapshot.expression_summary_default_cube.df[:]
    unique_genes = list(expr_default_cube.gene_ontology_term_id.unique())

    cube_query_params = WmgCubeQueryParams(
            cube_query_valid_attrs=READER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
            cube_query_valid_dims=READER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
        )
    
    query_obj = WmgQuery(snapshot, cube_query_params)
    profile_func(query_obj, unique_genes, randomize_reqs, num_trials, wait)

def print_profile_results_raw_query(snapshot, profile_func, randomize_reqs, num_trials, wait):
    expr_default_cube = snapshot.expression_summary_default_cube.df[:]
    unique_genes = list(expr_default_cube.gene_ontology_term_id.unique())

    cube_non_indexed_attrs = ['cell_type_ontology_term_id', 'disease_ontology_term_id', 'nnz', 'sum', 'sqsum']
    cube_indexed_dims = ['gene_ontology_term_id', 'tissue_ontology_term_id']
    profile_func(
        snapshot.expression_summary_cube,
        cube_indexed_dims,
        cube_non_indexed_attrs,
        unique_genes,
        randomize_reqs,
        num_trials,
        wait)

def gene_sample(randomize, all_genes, k):
    if randomize:
        return random.sample(all_genes, k=k)
    
    # If no randomization is required then just send a static list of 50 genes
    fifty_genes = ["ENSG00000188596","ENSG00000255794","ENSG00000130413",
                   "ENSG00000118997","ENSG00000165309","ENSG00000115423",
                   "ENSG00000197748","ENSG00000105877","ENSG00000186952",
                   "ENSG00000114670","ENSG00000186409","ENSG00000155761",
                   "ENSG00000152582","ENSG00000170959","ENSG00000157423",
                   "ENSG00000168038","ENSG00000174776","ENSG00000107249",
                   "ENSG00000253438","ENSG00000206530","ENSG00000133640",
                   "ENSG00000080298","ENSG00000162814","ENSG00000250305",
                   "ENSG00000156042","ENSG00000164692","ENSG00000168542",
                   "ENSG00000011465","ENSG00000108821","ENSG00000106809",
                   "ENSG00000139329","ENSG00000129009","ENSG00000113140",
                   "ENSG00000105894","ENSG00000248527","ENSG00000114115",
                   "ENSG00000100097","ENSG00000152661","ENSG00000234741",
                   "ENSG00000245910","ENSG00000087303","ENSG00000237550",
                   "ENSG00000175061","ENSG00000116774","ENSG00000189058",
                   "ENSG00000038427","ENSG00000198856","ENSG00000185222",
                   "ENSG00000203875","ENSG00000118418"]
    
    return fifty_genes

def get_query_k_random_genes_group_by_disease(randomize, all_genes, k=50):
    genes = gene_sample(randomize, all_genes, k)
    
    req = {
        "compare":"disease",
        "filter":{
            "dataset_ids":[],
            "development_stage_ontology_term_ids":[],
            "disease_ontology_term_ids":[],
            "gene_ontology_term_ids": genes,
            "organism_ontology_term_id":"NCBITaxon:9606",
            "self_reported_ethnicity_ontology_term_ids":[],
            "sex_ontology_term_ids":[],
            "publication_citations":[]
        },
        "is_rollup":True}
    compare = 'disease_ontology_term_id'
    
    criteria = WmgQueryCriteriaV2(**req["filter"])
    return criteria, compare

# Simple wall clock profiling using time module
def simple_wallclock_web_app_query(
        query_obj,
        unique_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=50):
    print(f"#####Simple Wallclock Profiling WEB APP QUERY#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}, randomize_queries: {randomize_reqs}, wait_time: {wait}")

    total_query_time = 0
    query_times = []
    for _ in range(num_trials):
        criteria, compare = get_query_k_random_genes_group_by_disease(
            randomize_reqs,
            unique_genes,
            num_genes_per_trial)
        start = time.perf_counter()
        result_df = query_obj.expression_summary(criteria, compare_dimension=compare)
        end = time.perf_counter()
        query_times.append(end - start)
        total_query_time += (end - start)
        time.sleep(wait)

    print(f"total_query_time: {total_query_time}, avg_query_time: {total_query_time/num_trials}")
    print(f"query_times: {query_times}")

def simple_wallclock_raw_query(
        expr_summary_array,
        cube_indexed_dims,
        cube_non_indexed_attrs,
        all_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=50):
    print(f"#####Simple Wallclock Profiling RAW QUERY#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}, randomize_queries: {randomize_reqs}, wait_time: {wait}")

    total_query_time = 0
    query_times = []
    for _ in range(num_trials):
        random_genes = gene_sample(randomize_reqs, all_genes, num_genes_per_trial)
        indexed_dims_query = tuple([random_genes, [], 'NCBITaxon:9606'])
        start = time.perf_counter()

        df_indexer_obj = expr_summary_array.query(
            cond=None,
            return_incomplete=True,
            use_arrow=True,
            attrs=cube_non_indexed_attrs,
            dims=cube_indexed_dims,
        ).df[indexed_dims_query]

        result_df = pd.concat(df_indexer_obj)

        end = time.perf_counter()
        query_times.append(end - start)
        total_query_time += (end - start)
        time.sleep(wait)

    print(f"total_query_time: {total_query_time}, avg_query_time: {total_query_time/num_trials}")
    print(f"query_times: {query_times}")

# Profiling cpu and non-cpu time using cProfile
def not_cpu_time():
    times = os.times()
    return times.elapsed - (times.system + times.user)


def cprofile_profiler_web_app_query(
        query_obj,
        unique_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=50,
        non_cpu_time=False):
    profiler = Profile(not_cpu_time) if non_cpu_time else Profile()
    
    for i in range(num_trials):
        print(f"###### Trial: {i+1} ######\n\n")
        criteria, compare = get_query_k_random_genes_group_by_disease(
            randomize_reqs,
            unique_genes,
            num_genes_per_trial)
        
        profiler.runcall(query_obj.expression_summary, criteria, compare_dimension=compare)
        s = io.StringIO()
        stats = Stats(profiler, stream=s)
        stats.strip_dirs()
        # 'tottime' is time spent in a function A EXCLUDING time spent in the functions A calls
        stats.sort_stats("tottime")
        stats.print_stats()

        # NOTE: This prints only the first 50 lines of the profile info.
        # To print full profile info do: `print(s.getvalue())``
        print('\n'.join(s.getvalue().split('\n')[:50]))
        print("\n\n\n")
        time.sleep(wait)

def cprofile_wallclock_web_app_query(
        query_obj,
        unique_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=50):
    print(f"#####CProfile Full Process Time Profiling WEB APP QUERY#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}, randomize_queries: {randomize_reqs}, wait_time: {wait}")
    cprofile_profiler_web_app_query(
        query_obj,
        unique_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=num_genes_per_trial)

def cprofile_non_cpu_time_web_app_query(
        query_obj,
        unique_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=50):
    print(f"#####CProfile Non CPU Time Profiling WEB APP QUERY#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}, randomize_queries: {randomize_reqs}, wait_time: {wait}")
    cprofile_profiler_web_app_query(
        query_obj,
        unique_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=num_genes_per_trial,
        non_cpu_time=True)

def cprofile_profiler_raw_query(
        expr_summary_array,
        cube_indexed_dims,
        cube_non_indexed_attrs,
        all_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=50,
        non_cpu_time=False):
    profiler = Profile(not_cpu_time) if non_cpu_time else Profile()

    def tiledb_raw_query(indexed_dims_query):
        df_indexer_obj = expr_summary_array.query(
            cond=None,
            return_incomplete=True,
            use_arrow=True,
            attrs=cube_non_indexed_attrs,
            dims=cube_indexed_dims,
        ).df[indexed_dims_query]

        result_df = pd.concat(df_indexer_obj)
        return result_df
    
    for i in range(num_trials):
        print(f"###### Trial: {i+1} ######\n\n")
        random_genes = gene_sample(randomize_reqs, all_genes, num_genes_per_trial)
        indexed_dims_query = tuple([random_genes, [], 'NCBITaxon:9606'])
        
        profiler.runcall(tiledb_raw_query, indexed_dims_query)
        s = io.StringIO()
        stats = Stats(profiler, stream=s)
        stats.strip_dirs()
        # 'tottime' is time spent in a function A EXCLUDING time spent in the functions A calls
        stats.sort_stats("tottime")
        stats.print_stats()

        # NOTE: This prints only the first 50 lines of the profile info.
        # To print full profile info do: `print(s.getvalue())``
        print('\n'.join(s.getvalue().split('\n')[:50]))
        print("\n\n\n")
        time.sleep(wait)

def cprofile_wallclock_raw_query(
        expr_summary_array,
        cube_indexed_dims,
        cube_non_indexed_attrs,
        all_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=50):
    print(f"#####CProfile Full Process Time Profiling RAW QUERY#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}, randomize_queries: {randomize_reqs}, wait_time: {wait}")
    cprofile_profiler_raw_query(
        expr_summary_array,
        cube_indexed_dims,
        cube_non_indexed_attrs,
        all_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=num_genes_per_trial)

def cprofile_non_cpu_time_raw_query(
        expr_summary_array,
        cube_indexed_dims,
        cube_non_indexed_attrs,
        all_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=50):
    print(f"#####CProfile Non CPU Time Profiling RAW QUERY#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}, randomize_queries: {randomize_reqs}, wait_time: {wait}")
    cprofile_profiler_raw_query(
        expr_summary_array,
        cube_indexed_dims,
        cube_non_indexed_attrs,
        all_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=num_genes_per_trial,
        non_cpu_time=True)

# Get profile information directly from tiledb stats dump
def tiledb_native_profiling_stats_web_app_query(
        query_obj, 
        unique_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=50):
    print(f"#####Tiledb Stats WEB APP QUERY#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}, randomize_queries: {randomize_reqs}, wait_time: {wait}")
    
    for i in range(num_trials):
        print(f"###### Trial: {i+1} ######\n\n")
        criteria, compare = get_query_k_random_genes_group_by_disease(
            randomize_reqs,
            unique_genes,
            num_genes_per_trial)
        tiledb.stats_enable()
        query_obj.expression_summary(criteria, compare_dimension=compare)
        tiledb_stats_str = tiledb.stats_dump(print_out=False, json=True)
        tiledb.stats_disable()
        pprint.pprint(tiledb_stats_str)
        print("\n\n\n")
        time.sleep(wait)

def tiledb_native_profiling_stats_raw_query(
        expr_summary_array,
        cube_indexed_dims,
        cube_non_indexed_attrs,
        all_genes,
        randomize_reqs,
        num_trials,
        wait,
        num_genes_per_trial=50):
    print(f"#####Tiledb Stats RAW QUERY#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}, randomize_queries: {randomize_reqs}, wait_time: {wait}")
    
    for i in range(num_trials):
        print(f"###### Trial: {i+1} ######\n\n")
        random_genes = gene_sample(randomize_reqs, all_genes, num_genes_per_trial)
        indexed_dims_query = tuple([random_genes, [], 'NCBITaxon:9606'])
        tiledb.stats_enable()
        df_indexer_obj = expr_summary_array.query(
            cond=None,
            return_incomplete=True,
            use_arrow=True,
            attrs=cube_non_indexed_attrs,
            dims=cube_indexed_dims,
        ).df[indexed_dims_query]
        result_df = pd.concat(df_indexer_obj)
        tiledb_stats_str = tiledb.stats_dump(print_out=False, json=True)
        tiledb.stats_disable()
        pprint.pprint(tiledb_stats_str)
        print("\n\n\n")
        time.sleep(wait)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Profile expression_summary cube query')
    
    parser.add_argument(
        '--data-location',
        required=True,
        choices=['local', 's3'],
        help='Specify data location (local or s3)')
    
    parser.add_argument(
        '--profile-type',
        required=True,
        choices=['simple-wallclock', 'cprofile-wallclock', 'cprofile-non-cpu', 'tiledb-stats'],
        help='Specify profile type')
    
    parser.add_argument(
        '--profile-scope',
        required=True,
        choices=['raw-query', 'web-app-query'],
        help='"raw-query" executes just the portion of the function that calls tiledb')
    
    parser.add_argument(
        '--num-trials',
        type=int,
        default=1,
        help='Specify number of trials to run')
    
    parser.add_argument(
        '--random',
        type=bool,
        default=False,
        help='Specify if the request should be randomly generated or not')
    
    parser.add_argument(
        '--wait',
        type=float,
        default=0.1,
        help='Specify number of seconds to wait between trials')
    
    args = parser.parse_args()

    if args.data_location == "local":
        snapshot = local_disk_snapshot()
    elif args.data_location == "s3":
        snapshot = s3_snapshot()
    
    # Profile the function executed by the web-app.
    # This involves constructing the arguments needed to
    # call the tiledb query function.
    # THE EXECUTION OF THIS FUNCTION REPRESENTS THE REAL
    # LATENCY OBSERVED BY THE USER
    if args.profile_scope == 'web-app-query':
        if args.profile_type == "simple-wallclock":
            profile_func = simple_wallclock_web_app_query
        elif args.profile_type == "cprofile-wallclock":
            profile_func = cprofile_wallclock_web_app_query
        elif args.profile_type == "cprofile-non-cpu":
            profile_func = cprofile_non_cpu_time_web_app_query
        elif args.profile_type == "tiledb-stats":
            profile_func = tiledb_native_profiling_stats_web_app_query
        
        print_profile_results_web_app_query(snapshot, profile_func, args.random, args.num_trials, args.wait)
    
    # Profile just the tiledb query. THE LATENCY OBSERVED HERE
    # IS ENTIRELY ATTRIBUTED TILED!!!
    if args.profile_scope == 'raw-query':
        if args.profile_type == "simple-wallclock":
            profile_func = simple_wallclock_raw_query
        elif args.profile_type == "cprofile-wallclock":
            profile_func = cprofile_wallclock_raw_query
        elif args.profile_type == "cprofile-non-cpu":
            profile_func = cprofile_non_cpu_time_raw_query
        elif args.profile_type == "tiledb-stats":
            profile_func = tiledb_native_profiling_stats_raw_query
        
        print_profile_results_raw_query(snapshot, profile_func, args.random, args.num_trials, args.wait)
    





