"""SPIKE helper: build chDB (embedded ClickHouse) MergeTree tables from the exported cube Parquet.

    PYTHONPATH=. venv-311/bin/python scripts/wmg_build_chdb.py <parquet_dir> <chdb_out_dir>

Reuses the Parquet produced by wmg_export_cube_to_parquet.py. For each cube it creates a MergeTree
table whose `ORDER BY` puts the always-filtered dims first (organism, then the selective gene id),
so `WHERE gene IN (...)` binary-searches the sparse primary index and reads only matching granules --
the layout property that makes selective lookups result-size instead of a scan. Secondary categorical
filters (cell_type, dataset, ...) get bloom_filter data-skipping indexes. The ORDER BY is THE tunable
this spike is testing; if pruning doesn't hold on the real cube, chDB fails the gate like Lance did.
"""

import os
import sys
import time

import pyarrow.parquet as pq
from chdb import session

CUBES = [
    "expression_summary",
    "expression_summary_default",
    "cell_counts",
    "marker_genes",
    "expression_summary_diffexp",
    "expression_summary_diffexp_simple",
    "cell_counts_diffexp",
]
DB = "wmg"

# PK per cube: low-cardinality organism first (cheap prune), then the selective gene id (binary search
# on the IN-list), then tissue. cell_counts/marker_genes have no gene dim; the diffexp cubes are keyed
# on the integer group_id, so ORDER BY group_id makes `group_id IN (...)` a binary-search slice.
ORDER_BY = {
    "expression_summary": ["organism_ontology_term_id", "gene_ontology_term_id", "tissue_ontology_term_id"],
    "expression_summary_default": ["organism_ontology_term_id", "gene_ontology_term_id", "tissue_ontology_term_id"],
    "cell_counts": ["organism_ontology_term_id", "tissue_ontology_term_id"],
    "marker_genes": ["organism_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"],
    "expression_summary_diffexp": ["group_id"],
    "expression_summary_diffexp_simple": ["group_id"],
    "cell_counts_diffexp": ["organism_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"],
}

# Secondary filter columns the API filters on but that aren't in the PK -> bloom_filter skip index.
SKIP_INDEX_COLS = [
    "cell_type_ontology_term_id",
    "tissue_original_ontology_term_id",
    "dataset_id",
    "disease_ontology_term_id",
    "self_reported_ethnicity_ontology_term_id",
    "sex_ontology_term_id",
    "publication_citation",
]


def _ch_type(pa_type) -> str:
    t = str(pa_type)
    if t in ("string", "large_string") or t.startswith("dictionary"):
        return "String"
    if t == "float":
        return "Float32"
    if t == "double":
        return "Float64"
    if t == "uint64":
        return "UInt64"
    if t == "uint32":
        return "UInt32"
    if t == "int64":
        return "Int64"
    if t == "int32":
        return "Int32"
    raise ValueError(f"unmapped arrow type {t}")


# junk columns the export's reset_index() leaves behind; not part of the cube schema.
IGNORE_COLS = {"index"}


def build_cube(sess, name: str, pq_path: str) -> int:
    schema = pq.read_schema(pq_path)
    cols = [(f.name, _ch_type(f.type)) for f in schema if f.name not in IGNORE_COLS]
    col_names = [c for c, _ in cols]
    order_by = [c for c in ORDER_BY[name] if c in col_names]
    skip = [c for c in SKIP_INDEX_COLS if c in col_names and c not in order_by]

    # backtick every identifier -- `index` and friends collide with ClickHouse keywords otherwise.
    col_defs = ",\n  ".join(f"`{c}` {t}" for c, t in cols)
    idx_defs = "".join(
        f",\n  INDEX skip_{c} `{c}` TYPE bloom_filter GRANULARITY 1" for c in skip
    )
    sess.query(f"DROP TABLE IF EXISTS {DB}.{name}")
    sess.query(
        f"CREATE TABLE {DB}.{name} (\n  {col_defs}{idx_defs}\n) "
        f"ENGINE = MergeTree ORDER BY ({', '.join(f'`{c}`' for c in order_by)})"
    )
    pq_abs = os.path.abspath(pq_path)
    select = ", ".join(f"`{c}`" for c in col_names)
    sess.query(f"INSERT INTO {DB}.{name} SELECT {select} FROM file('{pq_abs}', 'Parquet')")
    n = sess.query(f"SELECT count() FROM {DB}.{name}", "DataFrame").iloc[0, 0]
    return int(n)


def main(parquet_dir: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    sess = session.Session(os.path.join(out_dir, "db"))
    sess.query(f"CREATE DATABASE IF NOT EXISTS {DB}")
    try:
        for name in CUBES:
            pq_path = os.path.join(parquet_dir, f"{name}.parquet")
            if not os.path.isfile(pq_path):
                print(f"skip {name} (no parquet)")
                continue
            t = time.perf_counter()
            n = build_cube(sess, name, pq_path)
            print(f"{name}: {n:,} rows -> {DB}.{name} (ORDER BY {ORDER_BY[name]}) in {time.perf_counter()-t:.1f}s")
    finally:
        sess.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
