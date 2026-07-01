"""SPIKE: TileDB vs Lance on a REAL WMG cube snapshot.

    PYTHONPATH=. venv-311/bin/python scripts/wmg_lance_spike_realcube.py <snapshot_dir> <lance_dir>

Both backends read from local disk (TileDB array opened with tile cache; Lance dataset with scalar
indexes built by wmg_build_lance.py). Latencies are medians over many runs so the OS page cache warms
both sides. Parity asserted on expression_summary + cell_counts. Mirrors wmg_duckdb_spike_realcube.py
so the numbers drop into the same table.
"""

import sys
import time
from types import SimpleNamespace

import lance
import tiledb

from backend.common.census_cube.data.criteria import CensusCubeQueryCriteria
from backend.common.census_cube.data.query import CensusCubeQuery, CensusCubeQueryParams
from backend.common.census_cube.data.query_lance import LanceCensusCubeQuery, LanceCube
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


def open_lance(lance_dir):
    ns = SimpleNamespace()
    for name, attr, schema in CUBES:
        ds = lance.dataset(f"{lance_dir}/{name}.lance")
        setattr(ns, attr, LanceCube(ds, schema))
    return ns


def main(snapshot_dir, lance_dir):
    params = CensusCubeQueryParams(
        cube_query_valid_attrs=READER_CENSUS_CUBE_CUBE_QUERY_VALID_ATTRIBUTES,
        cube_query_valid_dims=READER_CENSUS_CUBE_CUBE_QUERY_VALID_DIMENSIONS,
    )
    t = time.perf_counter()
    tdb_ns = open_tiledb(snapshot_dir)
    print(f"opened TileDB cubes in {time.perf_counter()-t:.1f}s")
    lnc_ns = open_lance(lance_dir)

    tdb = CensusCubeQuery(tdb_ns, params)
    lnc = LanceCensusCubeQuery(lnc_ns, params)

    # pull real filter values straight off the Lance datasets
    es = lnc_ns.expression_summary_cube.dataset
    cc = lnc_ns.cell_counts_cube.dataset
    organism = cc.to_table(columns=["organism_ontology_term_id"], limit=1).to_pandas().iloc[0, 0]

    def distinct(ds, col, n, where=None):
        df = ds.to_table(columns=[col], filter=where).to_pandas()
        return df[col].dropna().unique().tolist()[:n]

    org_f = f"organism_ontology_term_id = '{organism}'"
    genes = distinct(es, "gene_ontology_term_id", 4, org_f)
    tissues = distinct(cc, "tissue_ontology_term_id", 3, org_f)
    cell_types = distinct(cc, "cell_type_ontology_term_id", 5, org_f)
    print(f"organism={organism}  genes={len(genes)} tissues={len(tissues)} cell_types={len(cell_types)}")

    def crit(**kw):
        return CensusCubeQueryCriteria(organism_ontology_term_id=organism, gene_ontology_term_ids=genes, **kw)

    cases = {
        "default (genes only)": (crit(), True),
        "genes + tissues": (crit(tissue_ontology_term_ids=tissues), False),
        "genes + cell_types (secondary)": (crit(cell_type_ontology_term_ids=cell_types), False),
    }

    print(f"\n{'case':<34} {'parity':<7} {'rows':>7} {'tiledb ms':>11} {'lance ms':>11} {'lance vs tdb':>13}")
    print("-" * 90)
    for label, (c, use_default) in cases.items():
        es_fn_tdb = (
            (lambda c=c: tdb.expression_summary_default(c)) if use_default else (lambda c=c: tdb.expression_summary(c))
        )
        es_fn_lnc = (
            (lambda c=c: lnc.expression_summary_default(c)) if use_default else (lambda c=c: lnc.expression_summary(c))
        )

        es_tdb, es_lnc = es_fn_tdb(), es_fn_lnc()
        par = parity("expression_summary", es_tdb, es_lnc)
        rows = len(normalize(es_tdb))

        tdb_ms = timeit(es_fn_tdb, n=5)
        lnc_ms = timeit(es_fn_lnc, n=5)
        flag = "PASS" if par.startswith("PASS") else "FAIL"
        ratio = lnc_ms / tdb_ms
        verdict = f"{ratio:.1f}x slower" if ratio >= 1 else f"{1/ratio:.1f}x faster"
        print(f"{label:<34} {flag:<7} {rows:>7} {tdb_ms:>11.2f} {lnc_ms:>11.2f} {verdict:>13}")
        if not par.startswith("PASS"):
            print(f"   {par}")

    cc_par = parity(
        "cell_counts", tdb.cell_counts(cases["genes + tissues"][0]), lnc.cell_counts(cases["genes + tissues"][0])
    )
    print(f"\ncell_counts parity: {cc_par}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
