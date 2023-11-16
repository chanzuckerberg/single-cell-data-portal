import argparse
import time
import random
import tiledb
import pandas as pd
import json

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

def simple_wallclock(snapshot, query_obj, num_trials=10, num_genes_per_trial=50):
    expr_default_cube = snapshot.expression_summary_default_cube.df[:]
    unique_genes = list(expr_default_cube.gene_ontology_term_id.unique())

    total_query_time = 0
    
    for _ in range(num_trials):
        criteria, compare = get_query_k_random_genes_group_by_disease(unique_genes, num_genes_per_trial)
        start = time.perf_counter()
        query_obj.expression_summary(criteria, compare_dimension=compare)
        end = time.perf_counter()
        total_query_time += (end - start)

    print(f"num_trials: {num_trials}, num_genes_per_trial: {num_genes_per_trial}")
    print(f"total_query_time: {total_query_time}, avg_query_time: {total_query_time/num_trials}")

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

    print_profile_results(snapshot, profile_func)
    





