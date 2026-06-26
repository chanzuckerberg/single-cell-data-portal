"""SPIKE helper: load the Parquet cubes into a persistent DuckDB database with native tables
and indexes on the slice columns — DuckDB's real "resident database" mode, the fair analog of
TileDB's opened+tiled array (vs a read_parquet VIEW that re-scans every query).

    PYTHONPATH=. python scripts/wmg_build_duckdb_db.py <parquet_dir> <db_path>
"""

import os
import sys
import time

import duckdb

# table -> columns to index (the dims queries slice on)
INDEXES = {
    "expression_summary": ["gene_ontology_term_id", "tissue_ontology_term_id"],
    "expression_summary_default": ["gene_ontology_term_id", "tissue_ontology_term_id"],
    "cell_counts": ["tissue_ontology_term_id"],
    "marker_genes": ["tissue_ontology_term_id", "cell_type_ontology_term_id"],
}


def main(parquet_dir, db_path):
    if os.path.exists(db_path):
        os.remove(db_path)
    con = duckdb.connect(db_path)
    con.execute("SET memory_limit='6GB'")
    con.execute(f"SET temp_directory='{os.path.dirname(db_path)}/.duckdb_tmp'")
    for table, cols in INDEXES.items():
        src = f"{parquet_dir}/{table}.parquet"
        t = time.perf_counter()
        con.execute(f"CREATE TABLE {table} AS SELECT * FROM read_parquet('{src}')")
        for col in cols:
            con.execute(f"CREATE INDEX idx_{table}_{col} ON {table}({col})")
        n = con.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
        print(f"{table}: {n:,} rows + {len(cols)} indexes in {time.perf_counter()-t:.1f}s")
    con.close()
    print(f"db size: {os.path.getsize(db_path)/1024/1024:.0f} MB")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
