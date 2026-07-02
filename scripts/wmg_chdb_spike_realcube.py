"""SPIKE: TileDB vs chDB (embedded ClickHouse MergeTree) on a REAL WMG cube snapshot.

    PYTHONPATH=. venv-311/bin/python scripts/wmg_chdb_spike_realcube.py <snapshot_dir> <chdb_dir>

Both backends read from local disk (TileDB array opened with tile cache; chDB MergeTree built by
wmg_build_chdb.py). Latencies are medians over many runs so the OS page cache warms both sides.
Parity asserted on expression_summary + cell_counts. Mirrors wmg_lance_spike_realcube.py so the
numbers drop into the same table.

Also runs the per-predicate diagnostic that exposed Lance's structural gap (gene-only vs organism-only
vs AND-combinations), with EXPLAIN indexes=1 to confirm the primary key actually prunes granules
instead of scanning -- the whole question this spike exists to answer.
"""

import os
import sys
import time
from types import SimpleNamespace

import tiledb
from chdb import session

from backend.common.census_cube.data.criteria import CensusCubeQueryCriteria
from backend.common.census_cube.data.query import CensusCubeQuery, CensusCubeQueryParams
from backend.common.census_cube.data.query_chdb import DB, ChdbCensusCubeQuery, ChdbCube
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


def open_chdb(chdb_dir):
    sess = session.Session(os.path.join(chdb_dir, "db"))
    ns = SimpleNamespace()
    for name, attr, schema in CUBES:
        setattr(ns, attr, ChdbCube(sess, name, schema))
    return sess, ns


def diagnostic(sess, genes, tissues, organism):
    """Per-predicate: rows matched, latency, and granule pruning (EXPLAIN) -- the Lance-killer test."""
    gene_in = "(" + ", ".join(f"'{g}'" for g in genes) + ")"
    tissue_in = "(" + ", ".join(f"'{t}'" for t in tissues) + ")"
    preds = {
        "gene only": f"gene_ontology_term_id IN {gene_in}",
        "organism only": f"organism_ontology_term_id = '{organism}'",
        "gene + organism": f"gene_ontology_term_id IN {gene_in} AND organism_ontology_term_id = '{organism}'",
        "gene + tissue": f"gene_ontology_term_id IN {gene_in} AND tissue_ontology_term_id IN {tissue_in}",
        "gene + tissue + organism": (
            f"gene_ontology_term_id IN {gene_in} AND tissue_ontology_term_id IN {tissue_in} "
            f"AND organism_ontology_term_id = '{organism}'"
        ),
    }
    print(f"\n{'filter':<28} {'rows':>12} {'latency ms':>11} {'granules':>16}")
    print("-" * 72)
    for label, where in preds.items():
        sql = f"SELECT count() FROM {DB}.expression_summary WHERE {where}"
        rows = int(sess.query(sql, "DataFrame").iloc[0, 0])
        ms = timeit(lambda sql=sql: sess.query(sql, "DataFrame"), n=5)
        exp = sess.query(f"EXPLAIN indexes=1 SELECT * FROM {DB}.expression_summary WHERE {where}", "DataFrame")
        txt = "\n".join(exp.iloc[:, 0].astype(str))
        gran = next((ln.strip() for ln in txt.splitlines() if "Granules:" in ln), "n/a")
        print(f"{label:<28} {rows:>12,} {ms:>11.1f} {gran:>16}")


def main(snapshot_dir, chdb_dir):
    params = CensusCubeQueryParams(
        cube_query_valid_attrs=READER_CENSUS_CUBE_CUBE_QUERY_VALID_ATTRIBUTES,
        cube_query_valid_dims=READER_CENSUS_CUBE_CUBE_QUERY_VALID_DIMENSIONS,
    )
    t = time.perf_counter()
    tdb_ns = open_tiledb(snapshot_dir)
    print(f"opened TileDB cubes in {time.perf_counter()-t:.1f}s")
    sess, ch_ns = open_chdb(chdb_dir)

    tdb = CensusCubeQuery(tdb_ns, params)
    ch = ChdbCensusCubeQuery(ch_ns, params)

    # pull real filter values straight off the MergeTree tables (LIMIT short-circuits the scan)
    def distinct(table, col, n, where=None):
        w = f" WHERE {where}" if where else ""
        df = sess.query(f"SELECT DISTINCT {col} FROM {DB}.{table}{w} LIMIT {n}", "DataFrame")
        return df[col].dropna().tolist()[:n]

    organism = sess.query(f"SELECT organism_ontology_term_id FROM {DB}.cell_counts LIMIT 1", "DataFrame").iloc[0, 0]
    org_f = f"organism_ontology_term_id = '{organism}'"
    genes = distinct("expression_summary", "gene_ontology_term_id", 4, org_f)
    tissues = distinct("cell_counts", "tissue_ontology_term_id", 3, org_f)
    cell_types = distinct("cell_counts", "cell_type_ontology_term_id", 5, org_f)
    print(f"organism={organism}  genes={len(genes)} tissues={len(tissues)} cell_types={len(cell_types)}")

    def crit(**kw):
        return CensusCubeQueryCriteria(organism_ontology_term_id=organism, gene_ontology_term_ids=genes, **kw)

    cases = {
        "default (genes only)": (crit(), True),
        "genes + tissues": (crit(tissue_ontology_term_ids=tissues), False),
        "genes + cell_types (secondary)": (crit(cell_type_ontology_term_ids=cell_types), False),
    }

    print(f"\n{'case':<34} {'parity':<7} {'rows':>7} {'tiledb ms':>11} {'chdb ms':>11} {'chdb vs tdb':>13}")
    print("-" * 90)
    for label, (c, use_default) in cases.items():
        es_fn_tdb = (
            (lambda c=c: tdb.expression_summary_default(c)) if use_default else (lambda c=c: tdb.expression_summary(c))
        )
        es_fn_ch = (
            (lambda c=c: ch.expression_summary_default(c)) if use_default else (lambda c=c: ch.expression_summary(c))
        )

        es_tdb, es_ch = es_fn_tdb(), es_fn_ch()
        par = parity("expression_summary", es_tdb, es_ch)
        rows = len(normalize(es_tdb))

        tdb_ms = timeit(es_fn_tdb, n=5)
        ch_ms = timeit(es_fn_ch, n=5)
        flag = "PASS" if par.startswith("PASS") else "FAIL"
        ratio = ch_ms / tdb_ms
        verdict = f"{ratio:.1f}x slower" if ratio >= 1 else f"{1/ratio:.1f}x faster"
        print(f"{label:<34} {flag:<7} {rows:>7} {tdb_ms:>11.2f} {ch_ms:>11.2f} {verdict:>13}")
        if not par.startswith("PASS"):
            print(f"   {par}")

    cc_par = parity(
        "cell_counts", tdb.cell_counts(cases["genes + tissues"][0]), ch.cell_counts(cases["genes + tissues"][0])
    )
    print(f"\ncell_counts parity: {cc_par}")

    diagnostic(sess, genes, tissues, organism)
    sess.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
