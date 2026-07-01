"""SPIKE helper: build Lance datasets CLUSTERED on the cube's indexed dims, + scalar indexes.

    PYTHONPATH=. venv-311/bin/python scripts/wmg_build_lance_sorted.py <parquet_dir> <lance_out_dir>

The plain (insertion-order) Lance spike regressed because a scalar-index lookup yields matching
row-ids that are *scattered* across 351M rows -> a costly take. This variant sorts each cube on its
TileDB dimension order before writing, so a `gene`/`gene+tissue` match is contiguous -> the take gets
cheap. i.e. we replicate TileDB's dimension tiling inside Lance. DuckDB does the out-of-core sort and
streams sorted batches straight into Lance (no intermediate parquet).
"""

import os
import sys
import time

import duckdb
import lance

from scripts.wmg_build_lance import INDEX_PLAN

# sort order = the cube's TileDB indexed dims, gene-leading (gene is always in the query criteria).
SORT_COLS = {
    "expression_summary": ["gene_ontology_term_id", "tissue_ontology_term_id", "organism_ontology_term_id"],
    "expression_summary_default": ["gene_ontology_term_id", "tissue_ontology_term_id", "organism_ontology_term_id"],
    "cell_counts": ["tissue_ontology_term_id", "organism_ontology_term_id"],
    "marker_genes": ["tissue_ontology_term_id", "cell_type_ontology_term_id", "organism_ontology_term_id"],
}


def build_cube(con, pq_path: str, lance_uri: str, sort_cols) -> int:
    order = ", ".join(sort_cols)
    rel = con.sql(f"SELECT * FROM read_parquet('{pq_path}') ORDER BY {order}")
    reader = rel.fetch_record_batch(rows_per_batch=500_000)  # DuckDB spills the sort; batches stream out sorted
    ds = lance.write_dataset(reader, lance_uri, schema=reader.schema, mode="overwrite", max_rows_per_file=2_000_000)
    cols = set(reader.schema.names)
    for col, idx_type in INDEX_PLAN.items():
        if col in cols:
            ds.create_scalar_index(col, index_type=idx_type)
    return ds.count_rows()


def main(parquet_dir: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    con = duckdb.connect()
    con.execute("SET memory_limit='6GB'")  # force the 351M-row sort to spill instead of OOM
    con.execute(f"SET temp_directory='{out_dir}/_duckdb_tmp'")
    for name, sort_cols in SORT_COLS.items():
        pq_path = os.path.join(parquet_dir, f"{name}.parquet")
        if not os.path.isfile(pq_path):
            print(f"skip {name} (no parquet)")
            continue
        lance_uri = os.path.join(out_dir, f"{name}.lance")
        t = time.perf_counter()
        n = build_cube(con, pq_path, lance_uri, sort_cols)
        print(f"{name}: {n:,} rows -> {lance_uri} (sorted on {sort_cols}, indexed) in {time.perf_counter()-t:.1f}s")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
