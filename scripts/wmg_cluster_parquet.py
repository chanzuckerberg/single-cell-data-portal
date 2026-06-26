"""SPIKE helper: cluster the exported Parquet by the slice dimensions so DuckDB can prune
row groups (the analog of TileDB tiling on indexed dims). Sorts expression_summary +
expression_summary_default by gene, tissue; copies the small cubes through.

    PYTHONPATH=. python scripts/wmg_cluster_parquet.py <in_dir> <out_dir>
"""

import os
import shutil
import sys
import time

import duckdb

SORT = {
    "expression_summary": "gene_ontology_term_id, tissue_ontology_term_id, organism_ontology_term_id",
    "expression_summary_default": "gene_ontology_term_id, tissue_ontology_term_id, organism_ontology_term_id",
}
COPY_THROUGH = ["cell_counts", "marker_genes"]


def main(in_dir, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    con = duckdb.connect()
    con.execute("SET memory_limit='6GB'")
    con.execute(f"SET temp_directory='{out_dir}/.duckdb_tmp'")
    for name, order_by in SORT.items():
        src, dst = f"{in_dir}/{name}.parquet", f"{out_dir}/{name}.parquet"
        t = time.perf_counter()
        con.execute(
            f"COPY (SELECT * FROM read_parquet('{src}') ORDER BY {order_by}) "
            f"TO '{dst}' (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 1000000)"
        )
        print(f"sorted {name} in {time.perf_counter()-t:.1f}s -> {os.path.getsize(dst)/1024/1024:.0f} MB")
    for name in COPY_THROUGH:
        shutil.copy(f"{in_dir}/{name}.parquet", f"{out_dir}/{name}.parquet")
        print(f"copied {name}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
