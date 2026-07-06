# WMG storage-backend spike — summary

**Question:** can the single-cell stack move off TileDB? **Answer:** it depends on the dataset's
*geometry*, and for the WMG cube specifically we now know both whether an alternative works (yes, chDB)
and how it would be deployed (a clickhouse-server sidecar). **Standing verdict: keep TileDB today** —
nothing forces a change, TileDB meets the latency budget, and the cube's upstream source is TileDB
regardless. This is a de-risking spike, not a migration.

Branch: `spike/wmg-parquet-duckdb` (pushed to `chanzuckerberg/single-cell-data-portal`, no PR, not
merged). Detail docs: [findings](STORAGE_BACKEND_MIGRATION_FINDINGS.md) ·
[forced-exit rework](WMG_CUBE_TILEDB_EXIT_REWORK.md) ·
[CH-vs-TileDB architecture](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md) ·
[census corpus generation](CENSUS_CORPUS_GENERATION_SCOPE.md).

---

## TL;DR

- **Two datasets, two answers.** The Explorer per-dataset matrix (`.cxg`) and the Census corpus are
  per-cell **tensors** → **Zarr fits** (`.cxg`→Zarr proven in the `single-cell-explorer` repo, PR
  #1369, ran in rdev). The WMG cube is a sparse predicate-filtered **OLAP aggregate** → different
  problem.
- **For the cube, the winning property is the storage *layout*, not the engine.** TileDB is fast
  because its data is clustered/tiled on the indexed dims, so a selective lookup is a contiguous slice
  (result-size access), not a scan. Any replacement must reproduce that.
- **Candidates, benchmarked on the real 351M-row cube:** DuckDB regressed 35–150× (scan-first); Lance
  regressed 1.2–10.7× (row-id/bitmap, not clustered — ~35 ms floor); **chDB (embedded ClickHouse
  MergeTree) PASSED** — parity everywhere, faster than TileDB on 4 of 5 shapes.
- **Concurrency:** embedded chDB is the wrong *serving* engine (dir lock, gevent serializes 14.7×), but
  a **clickhouse-server sidecar** fixes it (gevent overlaps, one shared ~300 MB process) and the
  read-only S3 `web`-disk serving model is validated locally.
- **Net:** a cube exit is viable end-to-end *if forced* — build with chDB, serve with a CH sidecar over
  a read-only S3 disk — but there's no driver to do it now, and it wouldn't remove TileDB from the
  deployment anyway (below).

---

## The investigation, in order

**1. Match the format to the geometry.** Zarr is a chunked N-d array store — no query engine, no
indexes, dense-native. It fits per-cell tensors (`.cxg`, corpus) but is wrong for the cube's
predicate-filtered OLAP shape. So for the cube the columnar candidate was the natural target.

**2. The gate.** Parity + latency on the **real v5 staging snapshot** (`1760292291`,
`expression_summary` = 351,721,698 rows). Median latency for every query shape must be **≤ TileDB's**
(6–15 ms). This is what flipped DuckDB from "2× faster" on a toy fixture to "150× slower" on the real
cube — nothing ships on theory.

**3. What lost, and why it's instructive.**

| Candidate | Result | Root cause |
|---|---|---|
| DuckDB / Parquet | 35–150× slower | scan-first; ART index unused for `IN`-lists |
| Lance (sorted + unsorted) | 1.2–10.7× slower, ~35 ms floor | scalar index returns row-ids, data not clustered → scattered take |
| Zarr | n/a | no query engine at all |

The lesson: it's not "having an index," it's whether the **layout clusters the leading dims**. TileDB
intersects predicates in multi-dimensional *coordinate* space; the losers intersect in row-id/bitmap
space.

**4. chDB passed.** MergeTree `ORDER BY (organism, gene, tissue)` gives a sparse primary index; a
`gene IN (...)` lookup binary-searches granule marks and reads only matching 8192-row granules —
**pruned to 3 of 42,969 granules**, the same result-size access TileDB gets. Parity PASS on every shape
(main + diffexp).

| shape | rows | TileDB | chDB | vs TileDB |
|---|---|---|---|---|
| default — genes only (60M) | 1,430 | 7.9 ms | 4.9 ms | **1.6× faster** |
| genes + tissues (351M) | 90 | 13.9 ms | 3.9 ms | **3.6× faster** |
| genes + cell_types (351M) | 103 | 8.5 ms | 10.0 ms | 1.2× slower |
| diffexp — full group_ids | 39,372 | 19.6 ms | 10.2 ms | **1.9× faster** |
| diffexp — simple group_ids | 39,372 | 16.1 ms | 6.4 ms | **2.5× faster** |

The one slower shape (`genes + cell_types`) is `cell_type` being off the sort-key prefix, served by a
bloom skip-index — marginal and tunable. Storage: chDB ~17 GB vs TileDB ~7.6 GB (tunable with
`LowCardinality`).

**5. Concurrency changed the deployment shape, not the verdict.** The latency numbers were
single-process. Prod is 5 gunicorn/gevent workers in a 16 GB task. A concurrency spike found embedded
chDB can't serve that model:

- exclusive **lock on its data dir** → workers can't share one on-disk DB;
- per-worker **copies** are wasteful (N×17 GB) and fragile (cloning a persisted dir fails to reopen);
- 10 gevent greenlets **serialized 14.7×** (a native call blocks the event loop).

**6. The serving fix: a clickhouse-server sidecar.** Build with embedded chDB → publish a read-only
MergeTree to S3 → serve with a `clickhouse-server` sidecar the workers hit over localhost. Validated
locally:

- **gevent overlaps** — 10 concurrent greenlets = **4.1×** (not 14.7×), because the query is now a
  socket call the hub can yield on;
- **one shared server at ~296 MB** serves all clients (vs 789 MB × N embedded); throughput 280–327 q/s
  at 5–10 clients (≥ TileDB's 5-worker 230 q/s);
- **read-only `web`-disk model works** — `static-files-disk-uploader` → serve over HTTP → attach
  read-only + LRU cache; cache **cold 52.6 ms → warm 4.5 ms**; two servers attach the *same* immutable
  data concurrently, no contention.
- Cost: per-query p50 ≈ 20 ms (HTTP transport overhead, not engine) — reducible via the native
  protocol (port 9000). Within the <10 s interactive budget.

**7. The residual — why this never fully removes TileDB.** The cube's *source* is the Census corpus,
read via `cellxgene_census.open_soma` → `tiledbsoma` (TileDB). The only TileDB-free artifact census
publishes is the per-dataset *pre-integration* source H5ADs; using them means re-implementing CZI's
integration (a fork of upstream ingestion). So even a perfect cube exit leaves TileDB in the pipeline
unless CZI republishes the integrated corpus in another format. Scoped in detail —
what the corpus build is, where its TileDB dependency actually sits (offline build-image, not
serving), and the options — in [`CENSUS_CORPUS_GENERATION_SCOPE.md`](CENSUS_CORPUS_GENERATION_SCOPE.md).

---

## Why chDB, not "a columnar DB"

The selection criterion is the storage **layout** (clustered + prefix pruning), not the "columnar"
label — which is exactly why two *columnar* candidates (DuckDB, Lance) lost. Engines that share the
pruning property (StarRocks, Doris, Pinot, Druid) are distributed clusters — heavier ops than the
"sync a dir from S3, no server" model. chDB/ClickHouse is distinctive in spanning embedded (build),
single-node sidecar (serve), and read-only-S3 web disk — it fits the existing deployment shape. The
one real semantic loss vs TileDB is the **single sort order** (one favored access path vs TileDB's
any-dimension coordinate tiling), mitigated by per-pattern tables / projections — which is already how
WMG organizes its 7 cubes. See the [architecture doc](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md).

---

## Recommendation

**Keep TileDB now.** No driver forces a change, the cube meets its budget, and a cube exit wouldn't
remove TileDB from the deployment (source read is upstream SOMA). The spike's value is **optionality**:
if a driver lands (dropping the C++ dependency, an object-store/ML-native mandate), the path is known
and de-risked, not a research project —

> **build with embedded chDB → publish a read-only MergeTree to S3 → serve with a clickhouse-server
> sidecar over a read-only `web`/`s3_plain` disk (+ local cache), workers as socket clients.**

**What's left before an exit is implementation, not research:** the 7 pipeline writers + snapshot
wiring; swap `_open_cube`/`create_ctx` for a CH client; then verify the two things a laptop can't —
**real-S3 cold-read latency** and behavior on real ECS hardware — plus switch to the native protocol
and tune the `cell_type` shape / `LowCardinality`. All bounded by the single `CensusCubeQuery._query`
seam (the DuckDB/Lance/chDB spikes are the proof the swap is contained).

---

## Artifacts

**Docs (repo root):** `STORAGE_BACKEND_MIGRATION_FINDINGS.md`, `WMG_CUBE_TILEDB_EXIT_REWORK.md`,
`CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md`, `CENSUS_CORPUS_GENERATION_SCOPE.md`, this summary.

**Backend seam (throwaway spikes):** `backend/common/census_cube/data/query_{chdb,duckdb,lance}.py`.

**Scripts (`scripts/`):** `wmg_export_cube_to_parquet.py`; chDB — `wmg_build_chdb.py`,
`wmg_chdb_spike_realcube.py`, `wmg_chdb_spike_diffexp.py`, `wmg_chdb_spike_concurrency.py`,
`wmg_ch_sidecar_concurrency.py`; DuckDB/Lance equivalents (`wmg_build_*`, `wmg_*_spike_realcube.py`,
`wmg_cluster_parquet.py`).

**Repro needs:** the spike venv (`venv-311/bin/python`), a local snapshot (interactive Okta login for
S3 — see the branch), `pip install chdb clickhouse-connect`, single-binary ClickHouse
(`curl https://clickhouse.com/ | sh`). Explorer `.cxg`→Zarr lives in the separate `single-cell-explorer`
repo (branch `poc/zarr-adaptor`, PR #1369).
