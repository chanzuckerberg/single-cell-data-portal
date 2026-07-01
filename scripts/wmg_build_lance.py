"""SPIKE helper: build Lance datasets (+ scalar indexes) from the exported cube Parquet.

    PYTHONPATH=. venv-311/bin/python scripts/wmg_build_lance.py <parquet_dir> <lance_out_dir>

Reuses the Parquet produced by wmg_export_cube_to_parquet.py. For each cube it writes a Lance
dataset, then builds the scalar indexes that make selective `col IN (...)` filters prune to
result-size: BTREE on the high-cardinality gene id, BITMAP on the low-cardinality categorical dims.
Without these indexes Lance scans like DuckDB did, so the index build is the whole point.
"""

import os
import sys
import time

import lance
import pyarrow.dataset as pads

CUBES = ["expression_summary", "expression_summary_default", "cell_counts", "marker_genes"]

# index_type per column when present. gene is ~30k distinct -> BTREE; the rest are low-cardinality
# categorical dims/attrs the API filters on -> BITMAP (roaring, result-size scaling).
INDEX_PLAN = {
    "gene_ontology_term_id": "BTREE",
    "organism_ontology_term_id": "BITMAP",
    "tissue_ontology_term_id": "BITMAP",
    "cell_type_ontology_term_id": "BITMAP",
    "tissue_original_ontology_term_id": "BITMAP",
    "dataset_id": "BITMAP",
    "disease_ontology_term_id": "BITMAP",
    "self_reported_ethnicity_ontology_term_id": "BITMAP",
    "sex_ontology_term_id": "BITMAP",
    "publication_citation": "BITMAP",
}


def build_cube(pq_path: str, lance_uri: str) -> int:
    # Stream Parquet -> Lance via a batch reader so the 4GB expression_summary never fully
    # materializes (mirrors the memory-bounded TileDB->Parquet export).
    pq_ds = pads.dataset(pq_path, format="parquet")
    reader = pq_ds.scanner(batch_size=500_000).to_reader()
    ds = lance.write_dataset(reader, lance_uri, schema=pq_ds.schema, mode="overwrite", max_rows_per_file=2_000_000)
    cols = set(pq_ds.schema.names)
    for col, idx_type in INDEX_PLAN.items():
        if col in cols:
            ds.create_scalar_index(col, index_type=idx_type)
    return ds.count_rows()


def main(parquet_dir: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    for name in CUBES:
        pq_path = os.path.join(parquet_dir, f"{name}.parquet")
        if not os.path.isfile(pq_path):
            print(f"skip {name} (no parquet)")
            continue
        lance_uri = os.path.join(out_dir, f"{name}.lance")
        t = time.perf_counter()
        n = build_cube(pq_path, lance_uri)
        mb = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(lance_uri) for f in fs) / 1024 / 1024
        print(f"{name}: {n:,} rows -> {lance_uri} ({mb:.0f} MB, indexed) in {time.perf_counter()-t:.1f}s")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
