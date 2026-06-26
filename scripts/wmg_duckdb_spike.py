"""SPIKE: compare TileDB vs Parquet+DuckDB for the WMG expression cube.

Runs representative /query-style criteria through BOTH backends behind the CensusCubeQuery seam,
asserts numeric parity, and times each. Uses the committed realistic-test-snapshot fixture.

    PYTHONPATH=. python scripts/wmg_duckdb_spike.py

Caveat: the fixture cube is tiny, so absolute latencies are not production numbers — the useful
signals here are (a) parity and (b) relative behaviour on identical data.
"""

import statistics
import tempfile
import time

import pandas as pd

from backend.common.census_cube.data.criteria import CensusCubeQueryCriteria
from backend.common.census_cube.data.query import CensusCubeQuery, CensusCubeQueryParams
from backend.common.census_cube.data.query_duckdb import DuckDBCensusCubeQuery, build_duckdb_snapshot
from backend.wmg.api.config import (
    READER_CENSUS_CUBE_CUBE_QUERY_VALID_ATTRIBUTES,
    READER_CENSUS_CUBE_CUBE_QUERY_VALID_DIMENSIONS,
)
from tests.unit.backend.wmg.fixtures.test_snapshot import FIXTURES_ROOT, load_realistic_test_snapshot

TEST_SNAPSHOT = "realistic-test-snapshot"

# Representative criteria. organism is required; the first hits the "default" (no-secondary-filter)
# cube path, the rest exercise indexed-dim slicing and non-indexed predicate filtering.
HUMAN = "NCBITaxon:9606"


def _criteria(**kw):
    return CensusCubeQueryCriteria(organism_ontology_term_id=HUMAN, **kw)


def _pick(snapshot, col, n=3):
    """Grab a few real values for a column from the fixture cell_counts to build valid filters."""
    return snapshot.cell_counts_df[col].dropna().unique().tolist()[:n]


def normalize(df):
    df = df.reset_index()
    df = df[[c for c in df.columns if c != "index"]]
    # TileDB ascii cols come back as bytes; DuckDB returns str. Decode so the two compare equal.
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].map(lambda v: v.decode() if isinstance(v, (bytes, bytearray)) else v)
    df = df.sort_values(list(df.columns)).reset_index(drop=True)
    for c in df.select_dtypes("float").columns:
        df[c] = df[c].round(4)
    return df


def parity(name, a, b):
    na, nb = normalize(a), normalize(b)
    cols = [c for c in na.columns if c in nb.columns]
    na, nb = na[cols], nb[cols]
    try:
        pd.testing.assert_frame_equal(na, nb, check_dtype=False, check_like=True, rtol=1e-4, atol=1e-6)
        return f"PASS ({len(na)} rows, cols={cols})"
    except AssertionError as e:
        return f"FAIL ({len(na)} vs {len(nb)} rows)\n    {str(e).splitlines()[0]}"


def timeit(fn, n=25):
    times = []
    for _ in range(n):
        t = time.perf_counter()
        fn()
        times.append((time.perf_counter() - t) * 1000)
    return statistics.median(times)


def main():
    params = CensusCubeQueryParams(
        cube_query_valid_attrs=READER_CENSUS_CUBE_CUBE_QUERY_VALID_ATTRIBUTES,
        cube_query_valid_dims=READER_CENSUS_CUBE_CUBE_QUERY_VALID_DIMENSIONS,
    )

    with load_realistic_test_snapshot(TEST_SNAPSHOT) as tdb_snapshot, tempfile.TemporaryDirectory() as pq_dir:
        t = time.perf_counter()
        duck_snapshot = build_duckdb_snapshot(TEST_SNAPSHOT, pq_dir, FIXTURES_ROOT)
        build_ms = (time.perf_counter() - t) * 1000

        tdb = CensusCubeQuery(tdb_snapshot, params)
        duck = DuckDBCensusCubeQuery(duck_snapshot, params)

        genes = (
            _pick(tdb_snapshot, "gene_ontology_term_id", 5)
            if "gene_ontology_term_id" in tdb_snapshot.cell_counts_df
            else []
        )
        # cell_counts has no gene dim; pull genes from expression_summary instead
        genes = (
            duck_snapshot.expression_summary_cube.con.execute(
                "SELECT DISTINCT gene_ontology_term_id FROM expression_summary LIMIT 5"
            )
            .df()["gene_ontology_term_id"]
            .tolist()
        )
        tissues = _pick(tdb_snapshot, "tissue_ontology_term_id", 2)
        cell_types = _pick(tdb_snapshot, "cell_type_ontology_term_id", 3)

        cases = {
            "default (genes only)": (_criteria(gene_ontology_term_ids=genes), True),
            "genes + tissues": (_criteria(gene_ontology_term_ids=genes, tissue_ontology_term_ids=tissues), False),
            "genes + cell_types (secondary)": (
                _criteria(gene_ontology_term_ids=genes, cell_type_ontology_term_ids=cell_types),
                False,
            ),
        }

        print(f"\nDuckDB store build (CSV->Parquet->views): {build_ms:.0f} ms\n")
        print(f"{'case':<34} {'parity':<8} {'tiledb ms':>10} {'duckdb ms':>10}")
        print("-" * 70)
        for label, (crit, use_default) in cases.items():
            es_tdb = tdb.expression_summary_default(crit) if use_default else tdb.expression_summary(crit)
            es_duck = duck.expression_summary_default(crit) if use_default else duck.expression_summary(crit)
            cc_tdb = tdb.cell_counts(crit)
            cc_duck = duck.cell_counts(crit)

            es_par = parity("expression_summary", es_tdb, es_duck)
            cc_par = parity("cell_counts", cc_tdb, cc_duck)

            es_tdb_ms = timeit(
                (lambda c=crit: tdb.expression_summary_default(c))
                if use_default
                else (lambda c=crit: tdb.expression_summary(c))
            )
            es_duck_ms = timeit(
                (lambda c=crit: duck.expression_summary_default(c))
                if use_default
                else (lambda c=crit: duck.expression_summary(c))
            )
            print(f"{label:<34} {'P' if es_par.startswith('PASS') else 'F':<8} {es_tdb_ms:>10.2f} {es_duck_ms:>10.2f}")
            print(f"  expression_summary: {es_par}")
            print(f"  cell_counts:        {cc_par}")
        print()


if __name__ == "__main__":
    main()
