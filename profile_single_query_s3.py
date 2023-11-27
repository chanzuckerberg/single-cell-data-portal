import argparse
import time
import pprint
import pandas as pd
import tiledb

def tiledb_query_params():
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
    
    cube_indexed_dims = [
        'gene_ontology_term_id',
        'tissue_ontology_term_id']
    
    cube_non_indexed_attrs = [
        'cell_type_ontology_term_id',
        'disease_ontology_term_id',
        'nnz',
        'sum',
        'sqsum']
    
    indexed_dims_query = tuple([fifty_genes, [], 'NCBITaxon:9606'])

    return cube_indexed_dims, cube_non_indexed_attrs, indexed_dims_query

def open_local_disk_tiledb_array():
    return tiledb.open('staging-snapshot/expression_summary')

def open_s3_tiledb_array():
    return tiledb.open(
            's3://cellxgene-wmg-staging/snapshots/v3/1701021778/expression_summary',
            config=tiledb.Config(
                {"vfs.s3.region": "us-west-2", 'py.init_buffer_bytes': 128*1024**2}))

def profile_query_time(tiledb_array):
    index_dims, non_indexed_atts, query = tiledb_query_params()
    start = time.perf_counter()

    df_indexer_obj = tiledb_array.query(
            cond=None,
            return_incomplete=True,
            use_arrow=True,
            attrs=non_indexed_atts,
            dims=index_dims,
        ).df[query]
    result_df = pd.concat(df_indexer_obj)
    
    end = time.perf_counter()

    query_time = end - start
    print(f"Query Time: {query_time}")


def profile_query_tiledb_stats(tiledb_array):
    index_dims, non_indexed_atts, query = tiledb_query_params()
    tiledb.stats_enable()

    df_indexer_obj = tiledb_array.query(
            cond=None,
            return_incomplete=True,
            use_arrow=True,
            attrs=non_indexed_atts,
            dims=index_dims,
        ).df[query]
    result_df = pd.concat(df_indexer_obj)
    tiledb_stats_str = tiledb.stats_dump(print_out=False, json=True)
    tiledb.stats_disable()
    pprint.pprint(tiledb_stats_str)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Profile single query in expression_summary cube')
    
    parser.add_argument(
        '--data-location',
        required=True,
        choices=['local', 's3'],
        help='Specify data location (local or s3)')
    
    parser.add_argument(
        '--profile-type',
        required=True,
        choices=['simple-wallclock', 'tiledb-stats'],
        help='Specify profile type')
    
    args = parser.parse_args()

    if args.data_location == "local":
        tiledb_array = open_local_disk_tiledb_array()
    elif args.data_location == "s3":
        tiledb_array = open_s3_tiledb_array()
    
    if args.profile_type == 'simple-wallclock':
        profile_query_time(tiledb_array)
    elif args.profile_type == 'tiledb-stats':
        profile_query_tiledb_stats(tiledb_array)
    





