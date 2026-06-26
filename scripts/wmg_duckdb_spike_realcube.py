"""SPIKE: TileDB vs Parquet+DuckDB on a REAL WMG cube snapshot.

    PYTHONPATH=. python scripts/wmg_duckdb_spike_realcube.py <snapshot_dir> <parquet_dir>

Both backends read from local disk (TileDB array opened with tile cache; DuckDB scanning Parquet
with predicate+projection pushdown). Latencies are medians over many runs, so the OS page cache
warms both sides — apples-to-apples. Parity asserted on expression_summary + cell_counts.
"""

import sys
import time
from types import SimpleNamespace

import duckdb
import tiledb

from backend.common.census_cube.data.criteria import CensusCubeQueryCriteria
from backend.common.census_cube.data.query import CensusCubeQuery, CensusCubeQueryParams
from backend.common.census_cube.data.query_duckdb import DuckDBCensusCubeQuery, DuckDBCube
from backend.common.census_cube.data.schemas.cube_schema import cell_counts_schema, expression_summary_schema
from backend.common.census_cube.data.schemas.cube_schema_default import (
    expression_summary_schema as expression_summary_default_schema,
)
from backend.common.census_cube.data.schemas.marker_gene_cube_schema import marker_genes_schema
from backend.common.census_cube.data.tiledb import create_ctx
from backend.wmg.api.config import (
    READER_CENSUS_CUBE_CUBE_QUERY_VALID_ATTRIBUTES,
    READER_CENSUS_CUBE_CUBE_QUERY_VALID_DIMENSIONS,
)
from scripts.wmg_duckdb_spike import normalize, parity, timeit

CUBES = [
    ("expression_summary", "expression_summary_cube", expression_summary_schema),
    ("expression_summary_default", "expression_summary_default_cube", expression_summary_default_schema),
    ("cell_counts", "cell_counts_cube", cell_counts_schema),
    ("marker_genes", "marker_genes_cube", marker_genes_schema),
]


def open_tiledb(snapshot_dir):
    ctx = create_ctx()
    ns = SimpleNamespace()
    for name, attr, _ in CUBES:
        setattr(ns, attr, tiledb.open(f"{snapshot_dir}/{name}", ctx=ctx))
    ns.cell_counts_df = ns.cell_counts_cube.df[:].reset_index()
    return ns


def open_duckdb(source):
    """source: a dir of Parquet files (VIEW over read_parquet, re-scans each query) OR a
    persistent .duckdb file (native indexed tables, resident — the fair analog of TileDB)."""
    ns = SimpleNamespace()
    if source.endswith(".duckdb"):
        con = duckdb.connect(source, read_only=True)
        for name, attr, schema in CUBES:
            setattr(ns, attr, DuckDBCube(con, name, schema))
    else:
        con = duckdb.connect()
        for name, attr, schema in CUBES:
            con.execute(f"CREATE VIEW {name} AS SELECT * FROM read_parquet('{source}/{name}.parquet')")
            setattr(ns, attr, DuckDBCube(con, name, schema))
    ns._con = con
    return ns


def main(snapshot_dir, parquet_dir):
    params = CensusCubeQueryParams(
        cube_query_valid_attrs=READER_CENSUS_CUBE_CUBE_QUERY_VALID_ATTRIBUTES,
        cube_query_valid_dims=READER_CENSUS_CUBE_CUBE_QUERY_VALID_DIMENSIONS,
    )
    t = time.perf_counter()
    tdb_ns = open_tiledb(snapshot_dir)
    print(f"opened TileDB cubes in {time.perf_counter()-t:.1f}s")
    duck_ns = open_duckdb(parquet_dir)

    tdb = CensusCubeQuery(tdb_ns, params)
    duck = DuckDBCensusCubeQuery(duck_ns, params)

    con = duck_ns._con
    organism = con.execute("SELECT organism_ontology_term_id FROM cell_counts LIMIT 1").fetchone()[0]
    genes = [
        r[0]
        for r in con.execute(
            f"SELECT DISTINCT gene_ontology_term_id FROM expression_summary "
            f"WHERE organism_ontology_term_id = '{organism}' LIMIT 4"
        ).fetchall()
    ]
    tissues = [
        r[0]
        for r in con.execute(
            f"SELECT DISTINCT tissue_ontology_term_id FROM cell_counts WHERE organism_ontology_term_id='{organism}' LIMIT 3"
        ).fetchall()
    ]
    cell_types = [
        r[0]
        for r in con.execute(
            f"SELECT DISTINCT cell_type_ontology_term_id FROM cell_counts WHERE organism_ontology_term_id='{organism}' LIMIT 5"
        ).fetchall()
    ]
    print(f"organism={organism}  genes={len(genes)} tissues={len(tissues)} cell_types={len(cell_types)}")

    def crit(**kw):
        return CensusCubeQueryCriteria(organism_ontology_term_id=organism, gene_ontology_term_ids=genes, **kw)

    cases = {
        "default (genes only)": (crit(), True),
        "genes + tissues": (crit(tissue_ontology_term_ids=tissues), False),
        "genes + cell_types (secondary)": (crit(cell_type_ontology_term_ids=cell_types), False),
    }

    print(f"\n{'case':<34} {'parity':<7} {'rows':>7} {'tiledb ms':>11} {'duckdb ms':>11} {'duck vs tdb':>12}")
    print("-" * 89)
    for label, (c, use_default) in cases.items():
        es_fn_tdb = (
            (lambda c=c: tdb.expression_summary_default(c)) if use_default else (lambda c=c: tdb.expression_summary(c))
        )
        es_fn_duck = (
            (lambda c=c: duck.expression_summary_default(c))
            if use_default
            else (lambda c=c: duck.expression_summary(c))
        )

        es_tdb, es_duck = es_fn_tdb(), es_fn_duck()
        par = parity("expression_summary", es_tdb, es_duck)
        rows = len(normalize(es_tdb))

        tdb_ms = timeit(es_fn_tdb, n=5)
        duck_ms = timeit(es_fn_duck, n=5)
        flag = "PASS" if par.startswith("PASS") else "FAIL"
        ratio = duck_ms / tdb_ms
        verdict = f"{ratio:.0f}x slower" if ratio >= 1 else f"{1/ratio:.0f}x faster"
        print(f"{label:<34} {flag:<7} {rows:>7} {tdb_ms:>11.2f} {duck_ms:>11.2f} {verdict:>12}")
        if not par.startswith("PASS"):
            print(f"   {par}")

    # cell_counts parity (no gene dim)
    cc_par = parity(
        "cell_counts", tdb.cell_counts(cases["genes + tissues"][0]), duck.cell_counts(cases["genes + tissues"][0])
    )
    print(f"\ncell_counts parity: {cc_par}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
