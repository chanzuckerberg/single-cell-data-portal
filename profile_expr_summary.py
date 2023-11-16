import argparse
import io
import os
import time
import random
import json
import pprint
from cProfile import Profile
from pstats import Stats

import tiledb
import pandas as pd

from backend.wmg.api.v2 import *

def local_disk_snapshot():
    return WmgSnapshot(snapshot_identifier="",
    expression_summary_cube=tiledb.open('staging-snapshot/expression_summary'),
    marker_genes_cube=tiledb.open('staging-snapshot/marker_genes'),
    expression_summary_default_cube=tiledb.open('staging-snapshot/expression_summary_default'),
    cell_counts_cube=tiledb.open('staging-snapshot/cell_counts'),
    cell_type_orderings=pd.read_json('staging-snapshot/cell_type_orderings.json').set_index(["tissue_ontology_term_id", "cell_type_ontology_term_id"])[
            "order"
        ].to_dict(),
    primary_filter_dimensions=json.load(open('staging-snapshot/primary_filter_dimensions.json','r')),
    dataset_metadata=json.load(open('staging-snapshot/dataset_metadata.json','r')), 
    filter_relationships=json.load(open('staging-snapshot/filter_relationships.json','r')))

def print_profile_results(snapshot, profile_func):
    cube_query_params = WmgCubeQueryParams(
            cube_query_valid_attrs=READER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
            cube_query_valid_dims=READER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
        )
    
    query_obj = WmgQuery(snapshot, cube_query_params)
    profile_func(snapshot, query_obj)

def get_query_k_random_genes_group_by_disease(all_genes, k=50):
    random_genes = random.sample(all_genes, k=k)
    
    req = {
        "compare":"disease",
        "filter":{
            "dataset_ids":[],
            "development_stage_ontology_term_ids":[],
            "disease_ontology_term_ids":[],
            "gene_ontology_term_ids": random_genes,
            "organism_ontology_term_id":"NCBITaxon:9606",
            "self_reported_ethnicity_ontology_term_ids":[],
            "sex_ontology_term_ids":[],
            "publication_citations":[]
        },
        "is_rollup":True}
    
    compare = req.get("compare", None)
    
    compare = find_dimension_id_from_compare(compare)
    
    criteria = WmgQueryCriteriaV2(**req["filter"])
    return criteria, compare

# Simple wall clock profiling using time module
def simple_wallclock(snapshot, query_obj, num_trials=10, num_genes_per_trial=50):
    print(f"#####Simple Wallclock Profiling#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}")

    expr_default_cube = snapshot.expression_summary_default_cube.df[:]
    unique_genes = list(expr_default_cube.gene_ontology_term_id.unique())

    total_query_time = 0
    for _ in range(num_trials):
        criteria, compare = get_query_k_random_genes_group_by_disease(unique_genes, num_genes_per_trial)
        start = time.perf_counter()
        query_obj.expression_summary(criteria, compare_dimension=compare)
        end = time.perf_counter()
        total_query_time += (end - start)

    print(f"total_query_time: {total_query_time}, avg_query_time: {total_query_time/num_trials}")

# Profiling cpu and non-cpu time using cProfile
def not_cpu_time():
    times = os.times()
    return times.elapsed - (times.system + times.user)


def cprofile_profiler(snapshot, query_obj, non_cpu_time=False, num_trials = 1, num_genes_per_trial=50):
    expr_default_cube = snapshot.expression_summary_default_cube.df[:]
    unique_genes = list(expr_default_cube.gene_ontology_term_id.unique())

    profiler = Profile(not_cpu_time) if non_cpu_time else Profile()
    for i in range(num_trials):
        print(f"###### Trial: {i+1} ######\n\n")
        criteria, compare = get_query_k_random_genes_group_by_disease(unique_genes, num_genes_per_trial)
        
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

def cprofile_wallclock(snapshot, query_obj, num_trials=1, num_genes_per_trial=50):
    print(f"#####CProfile Full Process Time Profiling#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}")
    cprofile_profiler(snapshot, query_obj, num_trials=num_trials, num_genes_per_trial=num_genes_per_trial)

def cprofile_non_cpu_time(snapshot, query_obj, num_trials=1, num_genes_per_trial=50):
    print(f"#####CProfile Non CPU Time Profiling#####")
    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}")
    cprofile_profiler(snapshot, query_obj, non_cpu_time=True, num_trials=num_trials, num_genes_per_trial=num_genes_per_trial)

# Get profile information directly from tiledb stats dump
def tiledb_native_profiling_stats(snapshot, query_obj, num_trials=1, num_genes_per_trial=50):
    expr_default_cube = snapshot.expression_summary_default_cube.df[:]
    unique_genes = list(expr_default_cube.gene_ontology_term_id.unique())

    for i in range(num_trials):
        print(f"###### Trial: {i+1} ######\n\n")
        criteria, compare = get_query_k_random_genes_group_by_disease(unique_genes, num_genes_per_trial)
        tiledb.stats_enable()
        query_obj.expression_summary(criteria, compare_dimension=compare)
        tiledb_stats_str = tiledb.stats_dump(print_out=False, json=True)
        tiledb.stats_disable()
        pprint.pprint(tiledb_stats_str)
        print("\n\n\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Profile expression_summary cube query on tiledb')

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
    
    args = parser.parse_args()

    if args.data_location == "local":
        snapshot = local_disk_snapshot()
    
    if args.profile_type == "simple-wallclock":
        profile_func = simple_wallclock
    elif args.profile_type == "cprofile-wallclock":
        profile_func = cprofile_wallclock
    elif args.profile_type == "cprofile-non-cpu":
        profile_func = cprofile_non_cpu_time
    elif args.profile_type == "tiledb-stats":
        profile_func = tiledb_native_profiling_stats

    print_profile_results(snapshot, profile_func)
    





