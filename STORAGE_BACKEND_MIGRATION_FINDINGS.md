# Storage backend migration: TileDB → Zarr / Parquet — findings

Investigation into moving off TileDB for two different datasets in the single-cell stack:

1. **Explorer per-dataset expression matrix** (`.cxg`) → **Zarr** — ✅ works, deployed to rdev.
2. **WMG / census expression cube** → **Zarr** (analysis) and **Parquet+DuckDB** (spiked) — ❌ don't migrate.

The headline: **match the storage format to the data's geometry.** A dense per-dataset tensor and a
sparse predicate-filtered OLAP cube are different shapes, and they want different backends — even
though both happen to use TileDB today.

---

## 1. Explorer `.cxg` → Zarr — WORKS

**Repo:** `single-cell-explorer`, branch `poc/zarr-adaptor` (PR #1369).

Explorer datasets are a per-dataset expression matrix + obs/var + embeddings — essentially an
**AnnData object**, a dense(ish) tensor. Zarr is the native fit (it's the chunked N-d array format,
and `anndata` reads/writes it directly).

**What was built:**
- `ZarrDataset(Dataset)` implementing the existing adaptor ABC via `anndata.read_zarr` — the REST,
  fbs-encoding, and compute layers were untouched (clean seam).
- One-line backend dispatch by `.zarr` suffix in `matrix_loader.py`.
- A `.cxg → .zarr` converter that round-trips through the existing adaptor.

**Result:** full parity vs the TileDB path (schema, obs/var, X slices, embeddings, diffexp — all
byte/value-identical), verified by unit tests **and end-to-end in rdev** (served a real
`pbmc3k.zarr` from S3 through the live Explorer UI).

**Friction (all surmountable):**
- Pinned `fsspec` 0.7.4 is too old for zarr's S3 `FSStore` array reads → the adaptor downloads the
  (small) zarr to a local temp dir and reads it locally. Upgrade path: bump fsspec/s3fs.
- pandas 2.x `future.infer_string` produces Arrow strings the zarr-v2 writer rejects → disabled in
  the converter.

**Verdict:** Zarr is a viable Explorer backend. Reads work; writes (user annotations) were out of
scope for the POC.

---

## 2. WMG / census cube — DON'T migrate

**Repo:** `single-cell-data-portal`, branch `spike/wmg-parquet-duckdb`.

The WMG cube is **not** a tensor. It's a **sparse, multi-dimensional, predicate-filtered OLAP
aggregate** (`expression_summary` ≈ 352M rows: gene × tissue × organism indexed dims, 7 non-indexed
filter attributes, `sum`/`sqsum`/`nnz`). The runtime API does **selective point lookups** — a handful
of genes, optionally sliced by tissue/cell_type — returning small result sets, with a <10s latency
budget. All TileDB access is isolated behind one seam (`CensusCubeQuery._query`).

### 2a. Why Zarr won't work for the cube (analysis)

Zarr is a chunked **N-dimensional array** store. It has **no query engine, no indexes, no predicate
pushdown**, and is **dense-native**. The cube's runtime depends on exactly the things Zarr lacks:

| Cube needs | TileDB provides | Zarr |
|---|---|---|
| Slice by dimension coordinate (`gene IN [...]`) | Indexed-dimension tiling | ✗ no index — must scan |
| Filter non-indexed attrs (`cell_type IN [...]`) | `QueryCondition` pushdown | ✗ nothing |
| Sparse storage (cube is sparse) | Native sparse arrays | ✗ dense; needs COO/CSR sidecar |
| Categorical compression | Dictionary filter | ✗ byte compression only |
| Dataframe-shaped I/O | `from_pandas` / `.df[]` | ✗ no table concept |

Putting the cube in raw Zarr means **reimplementing a query engine** (scan + pandas filtering) or
loading the whole cube into memory. Zarr is the wrong target for a relational/OLAP shape. The natural
columnar candidate is Parquet — so we spiked that instead.

### 2b. Parquet + DuckDB — spiked on the real cube, loses on latency

Built `DuckDBCensusCubeQuery` behind the same seam (criteria → SQL `IN` filters; reused
`CensusCubeQueryParams` so column projection matches by construction). Benchmarked on the **real v5
staging snapshot** (`expression_summary` = 351,721,698 rows) across three Parquet layouts.

**Parity: PASS everywhere** (expression_summary default/indexed/secondary + cell_counts, numerically
identical to TileDB).

**Latency (median, real cube):**

| Backend config | default (60M-row cube) | genes+tissues (351M-row) | genes+cell_types |
|---|---|---|---|
| **TileDB** (opened, tiled) | **6–10 ms** | **9–14 ms** | **6–15 ms** |
| DuckDB — Parquet `read_parquet` view | 327 ms | 1948 ms | 1567 ms |
| DuckDB — gene-sorted Parquet | 316 ms | 1870 ms | 1480 ms |
| DuckDB — native table + indexes | 217 ms | 1300 ms | 995 ms |

**DuckDB is 35–150× slower**, even sorted and indexed, despite tiny result sets (1–6k rows).

**Why:** latency scales with **table size**, not result size (611 MB → ~220 ms, 4.5 GB → ~1300 ms) —
the signature of a **full scan**. DuckDB is a scan-first analytical engine: its ART indexes aren't
used for `IN`-list filters and the planner prefers parallel scans. TileDB's array is **tiled/sorted on
its indexed dims**, so a few-gene lookup is a direct slice. WMG's access pattern is selective
point-lookups — an indexed-access workload, which is TileDB's strength and DuckDB's weakness.

**Storage got worse too:** TileDB 4.2 GB → Parquet 5.0 GB → native DuckDB DB **17.8 GB** (with indexes).

**Verdict:** Migrating the WMG cube to Parquet+DuckDB would regress interactive latency by 1–2 orders
of magnitude. Keep TileDB — its dimension tiling is the feature, not incidental. (Parquet+DuckDB would
only win for large-fraction scans/aggregations; WMG is not that.)

> Note: the *source* of the cube is `cellxgene-census` via `tiledbsoma` (itself TileDB, owned
> upstream). A full TileDB exit isn't possible without forking census ingestion regardless.

---

## The unifying principle

- **Dense per-dataset tensor (Explorer X / AnnData)** → **Zarr**. Slicing a hyperslab is the job; no
  query engine needed.
- **Sparse predicate-filtered OLAP cube with selective lookups (WMG)** → **keep TileDB**. The indexed
  dimension tiling is exactly what makes the interactive API fast; neither Zarr nor DuckDB replaces it.

The early intuition ("the cube is a columnar table → Parquet+DuckDB") was plausible but **wrong**, and
the spike on real data is what surfaced it. A small in-memory fixture initially showed DuckDB ~2×
*faster* — a misleading signal that only flipped at production scale.

---

## Artifacts

**Explorer (`single-cell-explorer`, branch `poc/zarr-adaptor`, PR #1369):**
- `server/dataset/zarr_dataset.py`, `server/dataset/matrix_loader.py` (dispatch),
  `scripts/cxg_to_zarr.py`, `server/tests/unit/dataset/test_zarr_dataset.py`.

**WMG (`single-cell-data-portal`, branch `spike/wmg-parquet-duckdb`):**
- `backend/common/census_cube/data/query_duckdb.py` — DuckDB backend behind the seam.
- `scripts/wmg_duckdb_spike.py` — fixture parity + latency.
- `scripts/wmg_export_cube_to_parquet.py` — stream TileDB cube → Parquet.
- `scripts/wmg_cluster_parquet.py` — gene-sort Parquet for row-group pruning.
- `scripts/wmg_build_duckdb_db.py` — load into native indexed DuckDB DB.
- `scripts/wmg_duckdb_spike_realcube.py` — real-snapshot benchmark.

**Reproduce the WMG benchmark** (needs a local cube snapshot dir + the spike venv):
```
PYTHONPATH=. python scripts/wmg_export_cube_to_parquet.py <snapshot_dir> <pq_dir>
PYTHONPATH=. python scripts/wmg_build_duckdb_db.py <pq_dir> <db.duckdb>
PYTHONPATH=. python scripts/wmg_duckdb_spike_realcube.py <snapshot_dir> <db.duckdb>
```

_These are throwaway spikes; not production code._
