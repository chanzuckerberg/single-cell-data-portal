"""SPIKE: TileDB vs chDB on the DIFFEXP query path (the group_id-keyed cube shape).

    PYTHONPATH=. venv-311/bin/python scripts/wmg_chdb_spike_diffexp.py <snapshot_dir> <chdb_dir>

Closes the gap left by wmg_chdb_spike_realcube.py, which only covered the categorical-dim cubes. The
diffexp path (backend/de/api) is different: it filters the small cell_counts_diffexp table to a set of
integer group_ids, then slices expression_summary_diffexp by `group_id IN (...)`. MergeTree ORDER BY
(group_id) should binary-search those groups. Asserts parity on the (expression, cell_counts) tuple
and times both backends for use_simple in {False, True}.
"""

import sys
import time
from types import SimpleNamespace

import pandas as pd
import tiledb

from backend.common.census_cube.data.criteria import BaseQueryCriteria
from backend.common.census_cube.data.query import CensusCubeQuery
from backend.common.census_cube.data.query_chdb import ChdbCensusCubeQuery, ChdbCube
from backend.common.census_cube.data.schemas.cube_schema_diffexp import expression_summary_schema
from backend.common.census_cube.data.tiledb import create_ctx
from scripts.wmg_duckdb_spike import normalize, parity, timeit

DIFFEXP_CUBES = [
    ("expression_summary_diffexp", "expression_summary_diffexp_cube"),
    ("expression_summary_diffexp_simple", "expression_summary_diffexp_simple_cube"),
]


def open_tiledb(snapshot_dir):
    ctx = create_ctx()
    ns = SimpleNamespace()
    for name, attr in DIFFEXP_CUBES:
        setattr(ns, attr, tiledb.open(f"{snapshot_dir}/{name}", ctx=ctx))
    with tiledb.open(f"{snapshot_dir}/cell_counts_diffexp", ctx=ctx) as cc:
        ns.cell_counts_diffexp_df = cc.df[:].reset_index()
    return ns


def open_chdb(chdb_dir, snapshot_dir, parquet_dir):
    from chdb import session

    sess = session.Session(f"{chdb_dir}/db")
    ns = SimpleNamespace()
    for name, attr in DIFFEXP_CUBES:
        setattr(ns, attr, ChdbCube(sess, name, expression_summary_schema))
    ccdf = pd.read_parquet(f"{parquet_dir}/cell_counts_diffexp.parquet")
    ns.cell_counts_diffexp_df = ccdf.drop(columns=["index"], errors="ignore")
    return sess, ns


def main(snapshot_dir, chdb_dir, parquet_dir="/tmp/wmg_pq"):
    t = time.perf_counter()
    tdb_ns = open_tiledb(snapshot_dir)
    print(f"opened TileDB diffexp cubes in {time.perf_counter()-t:.1f}s")
    sess, ch_ns = open_chdb(chdb_dir, snapshot_dir, parquet_dir)

    tdb = CensusCubeQuery(tdb_ns)
    ch = ChdbCensusCubeQuery(ch_ns)

    # realistic criteria: organism + one tissue, to bound the group_id set (DE filters this way)
    df = ch_ns.cell_counts_diffexp_df
    organism = df["organism_ontology_term_id"].iloc[0]
    tissue = df[df["organism_ontology_term_id"] == organism]["tissue_ontology_term_id"].iloc[0]
    crit = BaseQueryCriteria(organism_ontology_term_id=organism, tissue_ontology_term_ids=[tissue])
    print(f"organism={organism} tissue={tissue}")

    print(f"\n{'case':<22} {'parity':<7} {'groups':>7} {'rows':>9} {'tiledb ms':>11} {'chdb ms':>11} {'chdb vs tdb':>13}")
    print("-" * 92)
    for use_simple in (False, True):
        fn_tdb = lambda: tdb.expression_summary_and_cell_counts_diffexp(crit, use_simple)  # noqa: E731
        fn_ch = lambda: ch.expression_summary_and_cell_counts_diffexp(crit, use_simple)  # noqa: E731

        es_tdb, cc_tdb = fn_tdb()
        es_ch, cc_ch = fn_ch()
        # TileDB returns group_id already as a column; reset_index leaves a spurious RangeIndex "index"
        # column that chDB doesn't have -> drop it so parity compares the same 4 columns.
        es_tdb = es_tdb.reset_index().drop(columns=["index"], errors="ignore")
        key = "group_id_simple" if use_simple else "group_id"
        n_groups = cc_ch[key].nunique()

        par_es = parity("expression_diffexp", es_tdb, es_ch)
        par_cc = parity("cell_counts_diffexp", cc_tdb, cc_ch)
        rows = len(normalize(es_tdb))
        flag = "PASS" if par_es.startswith("PASS") and par_cc.startswith("PASS") else "FAIL"

        tdb_ms = timeit(fn_tdb, n=5)
        ch_ms = timeit(fn_ch, n=5)
        ratio = ch_ms / tdb_ms
        verdict = f"{ratio:.1f}x slower" if ratio >= 1 else f"{1/ratio:.1f}x faster"
        label = "simple group_ids" if use_simple else "full group_ids"
        print(f"{label:<22} {flag:<7} {n_groups:>7} {rows:>9} {tdb_ms:>11.2f} {ch_ms:>11.2f} {verdict:>13}")
        if flag == "FAIL":
            print(f"   es: {par_es}\n   cc: {par_cc}")

    sess.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
