# Storage backend migration: TileDB → Zarr / Parquet / Lance — findings

> For a one-page synthesis of the whole spike, see [`STORAGE_SPIKE_SUMMARY.md`](STORAGE_SPIKE_SUMMARY.md).

Investigation into moving off TileDB for two different datasets in the single-cell stack:

1. **Explorer per-dataset expression matrix** (`.cxg`) → **Zarr** — ✅ works, deployed to rdev.
2. **WMG / census expression cube** → **Zarr** (analysis), **Parquet+DuckDB** (spiked), **Lance**
   (spiked, sorted + unsorted, real cube) — ❌ these candidates regress; **chDB** is the one that passes.

The headline: **match the storage format to the data's geometry.** A dense per-dataset tensor and a
sparse predicate-filtered OLAP cube are different shapes, and they want different backends — even
though both happen to use TileDB today.

Migration-options status for the cube (§4): most non-TileDB backends were either analysed out (Zarr)
or benchmarked on the real 351M-row cube and **failed the no-regression gate** — DuckDB (35–150×
slower), Lance unsorted (1.2–10.7×), Lance clustered on the TileDB dim order (unchanged). **One
candidate passed: chDB (embedded ClickHouse MergeTree)** — spiked 2026-07-02 on the real cube, parity
on every shape and *faster* than TileDB on the main and diffexp paths (see §4). So the analysis
result: among tested backends, **chDB-embedded is the one proven-viable cube replacement** — the first
backend whose storage layout (MergeTree sparse primary index, granule-pruned `IN`-list lookups)
reproduces what makes TileDB fast — while DuckDB and Lance regress. Whether and when to act on that is
a product decision; the technical path is de-risked. Full rework scope in
[`WMG_CUBE_TILEDB_EXIT_REWORK.md`](WMG_CUBE_TILEDB_EXIT_REWORK.md); the architectural *why* (coordinate
tiling vs sorted-columnar prefix index) in
[`CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md`](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md).

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

## 2. WMG / census cube — only one backend passes the gate

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
of magnitude — so DuckDB is not a viable replacement. TileDB's dimension tiling is the feature here,
not incidental. (Parquet+DuckDB would only win for large-fraction scans/aggregations; WMG is not that.)

> Note: the *source* of the cube is `cellxgene-census` via `tiledbsoma` (itself TileDB, owned
> upstream). The corpus is geometrically a Zarr candidate (§3 — same AnnData/tensor shape as the
> `.cxg`), but that migration is CZI's to make, not the data portal's. From this repo's side the
> read leaves TileDB only two ways: (a) CZI republishes the integrated corpus in a non-TileDB format
> upstream, or (b) fork census ingestion from the per-dataset source H5ADs (`download_source_h5ad`) —
> the raw pre-integration inputs, which means re-harmonizing the corpus yourself. Neither is a
> storage-format swap in this repo, and there is no non-TileDB path to the *integrated* corpus today.

---

## 3. The three artifacts, two axes

All three single-cell artifacts in play sit on the **same TileDB storage engine**, but they split
along two *independent* axes: **API layer** (the SOMA data model vs. raw TileDB core) and **data
geometry** (per-cell tensor vs. OLAP aggregate). It's the geometry — not the API layer — that decides
whether Zarr fits.

| | **Census corpus** | **WMG cube** | **Explorer `.cxg`** |
|---|---|---|---|
| Produced by | Census team (upstream) | `backend/wmg/pipeline/` | `backend/layers/processing/process_cxg.py` |
| Input | — (the atlas itself) | Census corpus (SOMA) | labeled `local.h5ad` (per dataset) |
| Library | **`tiledbsoma`** | **`tiledb` core** | **`tiledb` core** |
| API layer | SOMA data model | raw arrays | raw TileDB group |
| Geometry | per-cell **tensor** (AnnData) | **OLAP aggregate** | per-cell **tensor** (AnnData) |
| Structure | `Experiment`: obs/var/X/obsm | standalone arrays (`expression_summary`, `cell_counts`…) | group: `X`, `obs`, `var`, `uns`, `emb/`, metadata |
| Granularity | cell × gene nonzero (billions) | (gene, tissue, organism, …) → `sum`/`sqsum`/`nnz` | cell × gene nonzero (one dataset) |
| X storage | SOMA sparse `X` layer | n/a (it *is* the aggregate) | sparse if <25% dense else dense; dims `obs`/`var` uint32, tiles 256/2048, byteshuffle+zstd-7, attr float32 |
| Index / dims | integer `soma_joinid` | `gene`/`tissue`/`organism` ontology IDs | **positional** `obs`/`var` (uint32 0..N) |
| obs/var | SOMA DataFrames | filter attrs live in the cube | dense 1-col TileDB arrays, zstd-22, categoricals as uint32 codes |
| Scope | whole atlas, both organisms | whole atlas, aggregated | **one dataset** |
| Access pattern | `axis_query` cell/gene slices | predicate-filtered aggregate point-lookups | random cell/gene **X slicing** + full obs load |
| Query engine | TileDB `QueryCondition` pushdown (obs/var `value_filter`) | TileDB `QueryCondition` + indexed-dim slicing — **load-bearing** | **none** — positional slicing only; compute runs in-process |
| Destination | `s3://cellxgene-census-public-us-west-2/...` | `s3://{WMG_BUCKET}/snapshots/v5/...` | `s3://{cellxgene_bucket}/{dataset_version_id}.cxg/` |
| **Zarr migration** | candidate (AnnData-native) | ❌ not viable for the cube | ✅ **proven** (Explorer `poc/zarr-adaptor`) |

The `.cxg` is what shows the two axes are independent: the corpus is tensor-shaped *and* SOMA; the
`.cxg` is tensor-shaped *but raw TileDB core* (no SOMA model — just positional `obs`/`var` dims and a
float32 X attr). So "TileDB-SOMA vs plain TileDB" is **not** the same distinction as "tensor vs OLAP."
The `.cxg` is the lower-level cousin of the corpus: same AnnData geometry, per-dataset, hand-built on
core TileDB for Explorer's interactive matrix slicing — which is exactly why it migrated to Zarr cleanly.

> Note: the cube's *source* is the Census corpus, read via `cellxgene_census.open_soma(census_version="latest")`.
> Resolution is plain JSON-over-HTTPS (`census.cellxgene.cziscience.com/.../release.json`), but it returns a
> `soma` locator (`s3://cellxgene-data-public/cell-census/<date>/soma/`, us-west-2, CZI-owned) that opens through
> `tiledbsoma` — i.e. TileDB. `get_anndata` is *not* a non-TileDB path; it reads through `tiledbsoma` too. The one
> TileDB-free artifact census offers is the per-dataset source H5ADs (`.../h5ads/<dataset_id>.h5ad`), but those are
> pre-integration raw inputs. The corpus's tensor geometry means it *is* a Zarr candidate (the "candidate" cell
> below), but that migration is upstream (CZI's) — this repo can only reach a non-TileDB read by consuming a
> future CZI republish or by forking ingestion from the source H5ADs. Neither is a storage-format swap here.

> External corroboration: CZI's [`chanzuckerberg/multimodal-slicing`](https://github.com/chanzuckerberg/multimodal-slicing)
> spike (Nov 2025) independently validates the corpus→Zarr candidate cell. It converts ~662 census-derived
> per-dataset h5ads to Zarr on S3 and serves lazy cell/gene slicing with obs-metadata filtering — the same
> **tensor geometry and positional-slicing access pattern** as the `.cxg` (no query engine; in-process compute),
> confirming the tensor→Zarr thesis on a second dataset. It also builds a shared "gene universe" (79,948 genes)
> to remap each dataset's local gene indices into a common space — i.e. a partial reimplementation of census's
> upstream harmonization, which is exactly the "fork ingestion" route (b) above, no longer hypothetical. Note it
> is a small-scale, I/O-bound demo (second-scale queries) and **makes no claim on the OLAP cube's workload** —
> it corroborates the Zarr half of this doc, not the cube findings (where common columnar engines regress).

---

## 4. If we *did* migrate the cube off TileDB — options that won't regress

§2 found that among common columnar engines only chDB matches TileDB on the cube. If a cube exit is
pursued (e.g. to drop the C++ dependency, or to align the cube with an object-store/ML-native stack),
these are the only options that can avoid a latency regression. Investigated 2026-06-30.

> Full rework scope for a forced cube exit — layers, touchpoints, per-layer rework, benchmarks, and
> testing — lives in the sibling doc [`WMG_CUBE_TILEDB_EXIT_REWORK.md`](WMG_CUBE_TILEDB_EXIT_REWORK.md).

### The reframe: replace the *layout property*, not "the engine"

TileDB's 6–15 ms didn't come from a query planner — it came from the data being **sorted/tiled on the
indexed dims**, so a few-gene lookup is a direct slice and latency scales with *result* size, not
*table* size. DuckDB lost because it's scan-first (ART index unused for `IN`-lists). So the
non-negotiable requirement for any replacement is: **the storage layout itself must index the leading
dims** so the engine reads only matching rows. That filter eliminates most candidates up front:

| Ruled OUT | Why |
|---|---|
| DuckDB / Parquet-as-table | proven 35–150× slower — full scan (§2b) |
| chDB querying raw Parquet | same scan trap — chDB only wins on its own MergeTree storage, not on Parquet |
| Zarr | no query engine at all (§2a) |
| In-memory pandas scan | O(n) mask, scales with table size |

### Viable options

**A. ClickHouse MergeTree — embedded via `chDB`** ✅ *spiked on the real cube 2026-07-02; PASSED the gate — the only candidate that does*
- *Engine replacement:* SQL. `ORDER BY (organism, gene, tissue)` → sparse primary index; `WHERE gene
  IN (...)` binary-searches in-memory index marks and reads only matching 8192-row granules. Filter
  dims get bloom data-skipping indexes. `sum`/`sqsum`/`nnz` already precomputed → query is a projected
  filter. Diffexp cubes keyed `ORDER BY (group_id)`.
- *Perf (real 351M cube, parity PASS everywhere):*

  | shape | TileDB | chDB | vs TileDB |
  |---|---|---|---|
  | default — genes only (60M) | 7.9 ms | 4.9 ms | **1.6× faster** |
  | genes + tissues (351M) | 13.9 ms | 3.9 ms | **3.6× faster** |
  | genes + cell_types (secondary, 351M) | 8.5 ms | 10.0 ms | 1.2× slower |
  | diffexp — full group_ids | 19.6 ms | 10.2 ms | **1.9× faster** |
  | diffexp — simple group_ids | 16.1 ms | 6.4 ms | **2.5× faster** |

  Granule diagnostic: gene-selective lookups prune to **3 / 42,969 granules** (~2–3 ms) — the
  result-size access TileDB has and Lance's row-id/bitmap model lacked (Lance had a ~35 ms floor).
- *Scope:* `chDB` embedded keeps the in-process, dir-to-S3 snapshot model — no server, no terraform
  change. (A ClickHouse *server* would break the dir-to-S3 model; not needed.)
- *Caveats:* secondary `cell_type` filter is the one shape ~1.2× slower (marginal, tunable via a `set`
  index or adding cell_type to the sort key); DB is ~17 GB vs TileDB ~7.6 GB (~2.2×).
- *Concurrency:* a follow-up spike showed embedded chDB is the wrong *serving* engine under load
  (exclusive dir lock → can't share a dir across workers; in-process gevent queries serialize ~15×).
  The fix is to **build with chDB, serve with a clickhouse-server sidecar over a read-only `web`/
  `s3_plain` disk** — queries become socket calls (gevent-friendly), caches shared, backpressure
  enforced. Full analysis in the rework doc §9.

**B. Lance** (embedded, object-store-native) — *spiked on the real cube; closest candidate but still regresses (details below)*
- *Engine replacement:* Lance scanner + scalar indexes — `BTREE` on `gene`, `BITMAP` on low-cardinality
  dims (`organism`/`tissue`/`cell_type`). Lance's bitmap search time "scales linearly with the number of
  values the query requires" — result-size, not table-size. Exactly the property we need.
- *Perf:* should match/beat the point-lookup; aggregates aren't engine-precomputed but result sets are
  1–6k rows (trivial). Must be benchmarked — same "looks right in theory" spot DuckDB fooled us in.
- *Scope:* medium; keeps today's architecture (embedded, no server, dir-to-S3 snapshot), aligns with the
  scverse/FM-training direction (same ecosystem pull as the zarr ask).
- *Risk:* OLAP-aggregate maturity; verify scalar-index pruning empirically.

**Spiked on the real cube (2026-06-30) — better than DuckDB, still regresses, fails the gate.**
Built `LanceCensusCubeQuery` behind the seam (`query_lance.py`), exported the real v5 staging snapshot
(`1760292291`, `expression_summary` = 351,721,698 rows) TileDB → Parquet → Lance, with BTREE on
`gene_ontology_term_id` + BITMAP on the low-cardinality dims. **Parity: PASS everywhere.**

| case | parity | rows | TileDB | Lance | vs TileDB |
|---|---|---|---|---|---|
| default — genes only (60M cube) | PASS | 8,691 | 8.1 ms | 10.0 ms | 1.2× slower |
| genes + tissues (351M) | PASS | 2,072 | 11.8 ms | 70.8 ms | 6.0× slower |
| genes + cell_types (351M) | PASS | 3,122 | 11.5 ms | 123.4 ms | 10.7× slower |

Lance is a large improvement over DuckDB (35–150× → 1.2–10.7×) and hits near-parity on the small
default cube, but **regresses on the 351M-row cube → fails the ≤-TileDB gate.** Why: the scalar index
returns matching row-ids cheaply, but the data isn't *clustered* on those dims, so Lance does a
scattered "take" of thousands of rows spread across 351M, vs TileDB's **contiguous tile slice** on its
sorted `(gene, tissue, organism)` dims. Regression scales with match scatter (secondary `cell_type`
filter is worst). The lesson sharpens: **the sort/clustering layout — not merely having an index — is
what makes TileDB fast.**

**Sorted/clustered Lance (2026-06-30) — sorting did NOT help; the gap is structural.** Rebuilt the
cube as Lance clustered on the TileDB dim order `(gene, tissue, organism)` (DuckDB out-of-core sort
streamed into Lance) + the same scalar indexes, then re-benchmarked:

| case | rows | TileDB | unsorted Lance | **sorted Lance** |
|---|---|---|---|---|
| default — genes only (60M) | 8,691 | ~7 ms | 10.0 ms | 8.9 ms |
| genes + tissues (351M) | 2,072 | ~12 ms | 70.8 ms | 69.6 ms |
| genes + cell_types (351M) | 3,122 | ~11 ms | 123.4 ms | 121.5 ms |

Clustering moved nothing. A per-predicate diagnostic on the sorted 351M cube shows why:

| filter | rows matched | latency |
|---|---|---|
| gene only | 13,979 | 35.5 ms |
| organism only | 319,599,548 | 11,251 ms |
| gene + organism | 13,979 | 44.1 ms |
| gene + tissue | 396 | 37.6 ms |
| gene + tissue + organism | 396 | 42.3 ms |

Two structural facts, neither tunable by layout:
1. **~35–44 ms fixed floor** for any gene-selective query — does not scale down with result size
   (396 rows costs the same as 14k). TileDB does the same lookups in 6–12 ms → Lance is ~3–6× slower
   even at its best.
2. **Low-cardinality dims produce dataset-sized bitmaps.** `organism` alone matches 319M rows (11 s, a
   roaring bitmap over half the table); survivable only because `gene` is ANDed first, but the secondary
   `cell_type` filter pays a bitmap-intersection cost that scales with **table size, not result size**.

Root cause: TileDB intersects predicates in **multi-dimensional coordinate space** and addresses the
exact cells; Lance intersects in **row-id/bitmap space** (fixed floor + low-card-bitmap cost). That's an
architectural difference, not a knob. **Verdict: Lance can't meet the no-regression gate — not a
viable cube replacement.** Of all candidates Lance is the closest (1.2–10.7× vs DuckDB's 35–150×), but
on this access pattern TileDB's multi-dim tiling has no equivalent; chDB is the one that closes the gap.

Artifacts: `backend/common/census_cube/data/query_lance.py`, `scripts/wmg_build_lance.py`,
`scripts/wmg_build_lance_sorted.py`, `scripts/wmg_lance_spike_realcube.py`.

**C. DIY: partition by gene + manifest index** — rebuild the tiling yourself. Guaranteed result-size
access, but it's literally reimplementing TileDB's dimension tiling in app code — if you're rebuilding
TileDB, TileDB already does it better. Only justified if dropping the C++ dep is the *sole* goal — and
A/B do it with less code you own.

### Scope (common to any choice — the seam bounds it)

The DuckDB spike proved the swap is contained to one class. Five surfaces:

1. **Query** — reimplement `CensusCubeQuery._query`. Must cover both shapes: `expression_summary`
   (gene/tissue/organism) *and* the diffexp cubes' integer `group_id` key.
2. **Write** — 7 array writers in `backend/wmg/pipeline/` → new-format writers.
3. **Snapshot/load** — `snapshot.py` open + S3 sync + `latest_snapshot_identifier` versioning. Lance/DIY
   keep the dir-to-S3 model; ClickHouse-server breaks it (chDB-embedded avoids this — biggest scope delta).
4. **Gate** — reuse the real-cube parity + latency harness from the spike.
5. **Source unchanged** — still reads TileDB-SOMA census. Never a full TileDB exit, whatever you pick.

### The hard gate

Before committing to any candidate: run the existing harness on the 351M-row staging snapshot and
require median latency for all three query shapes to be **≤ TileDB's**. That benchmark is what flipped
DuckDB from "2× faster" (toy fixture) to "150× slower" (real cube) — nothing ships on theory.

---

## The unifying principle

- **Dense per-dataset tensor (Explorer X / AnnData)** → **Zarr**. Slicing a hyperslab is the job; no
  query engine needed.
- **Sparse predicate-filtered OLAP cube with selective lookups (WMG)** → **TileDB or an equivalent
  sorted-columnar store (chDB)**. The indexed dimension tiling is what makes the interactive API fast;
  Zarr and DuckDB don't replace it — only chDB's MergeTree reproduces it.

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
- `backend/common/census_cube/data/query_chdb.py` — **chDB (MergeTree) backend behind the seam (the one that passed)**.
- `scripts/wmg_build_chdb.py` — Parquet → MergeTree (ORDER BY primary index + bloom skip-indexes).
- `scripts/wmg_chdb_spike_realcube.py` — real-snapshot benchmark + granule-pruning diagnostic.
- `scripts/wmg_chdb_spike_diffexp.py` — diffexp-path parity + latency.

**Reproduce the WMG benchmark** (needs a local cube snapshot dir + the spike venv):
```
PYTHONPATH=. python scripts/wmg_export_cube_to_parquet.py <snapshot_dir> <pq_dir>
# DuckDB (regresses):
PYTHONPATH=. python scripts/wmg_build_duckdb_db.py <pq_dir> <db.duckdb>
PYTHONPATH=. python scripts/wmg_duckdb_spike_realcube.py <snapshot_dir> <db.duckdb>
# chDB (passes):
PYTHONPATH=. python scripts/wmg_build_chdb.py <pq_dir> <chdb_dir>
PYTHONPATH=. python scripts/wmg_chdb_spike_realcube.py <snapshot_dir> <chdb_dir>
PYTHONPATH=. python scripts/wmg_chdb_spike_diffexp.py <snapshot_dir> <chdb_dir>
```

_These are throwaway spikes; not production code._
