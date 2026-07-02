# Rebuilding the WMG cube query engine to exit TileDB — rework scope

Sibling to [`STORAGE_BACKEND_MIGRATION_FINDINGS.md`](STORAGE_BACKEND_MIGRATION_FINDINGS.md). That
doc's verdict stands: **keep TileDB for the WMG cube** — every spiked backend regresses. This doc
answers a *different* question: **if a TileDB exit for the cube were forced anyway, what would it
actually take** to rebuild the query engine and every supporting layer? It weighs all options and
inventories the layers, indexing, rework scope, benchmarks, testing, and touchpoints.

---

## 1. Context & scope

**Premise:** a future driver forces a TileDB exit for the cube (e.g. dropping the TileDB C++
dependency from the deployment, or an org-wide object-store/ML-native storage mandate). Nothing
forces this today — see the recommendation in §9.

**In scope — the cube's own layers:**
- **Query/read** — the `CensusCubeQuery` seam and its consumers.
- **Write/build** — the 7 TileDB array writers in the WMG pipeline.
- **Snapshot/load/deploy** — how the cube is opened, versioned, synced, and served.

**Out of scope — the upstream source read (residual, stays TileDB):** the pipeline reads the
Census corpus via `cellxgene_census.open_soma(...)` → `tiledbsoma`, which *is* TileDB. Removing
that isn't a storage-format swap in this repo — it means forking CZI's census integration (see
§8 and the findings doc's §3 footnotes). Every option below still reads TileDB-SOMA at ingest.
There is no full TileDB exit; there is a **cube** exit.

---

## 2. What makes TileDB fast (the thing a replacement must reproduce)

> Deeper architectural comparison (coordinate tiling vs sorted-columnar prefix index, sparsity, query
> engine, deployment): [`CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md`](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md).

The cube's 6–15 ms interactive latency does **not** come from a query planner. It comes from the
data being **sorted/tiled on its indexed dimensions** (`gene`, `tissue`, `organism`), so a
few-gene point-lookup is a **direct contiguous slice** and latency scales with *result* size, not
*table* size. The runtime access pattern is selective point-lookups returning 1–6k rows out of
352M, under a <10s budget.

The spikes proved this is a **layout property, not an "index" property**:
- DuckDB with ART indexes → 35–150× slower (planner prefers parallel scans; ART unused for
  `IN`-lists).
- Lance with BTREE+BITMAP scalar indexes → still 1.2–10.7× slower. The index returns matching
  row-ids cheaply, but the data isn't *clustered* on those dims, so Lance does a scattered "take"
  across 351M rows. Sorting/clustering Lance on the TileDB dim order moved **nothing** — a ~35 ms
  fixed floor per gene-selective query, and low-cardinality dims (`organism` alone matches 319M
  rows) produce dataset-sized bitmaps.

**Root cause:** TileDB intersects predicates in **multi-dimensional coordinate space** and
addresses the exact cells; the columnar candidates intersect in **row-id/bitmap space** (fixed
floor + low-cardinality-bitmap cost that scales with table size). That's architectural.

**Therefore the non-negotiable requirement for any replacement:** the *storage layout* must index
and cluster the leading dims so the engine reads only matching rows — result-size access, not a
scan. This filter eliminates most candidates up front (§4).

---

## 3. The three layers + touchpoint inventory

The exit is contained behind three seams. Paths/lines are current as of this writing.

### Layer 1 — Query / read seam (`backend/common/census_cube/data/`)

| Element | Location | Notes |
|---|---|---|
| Seam class | `query.py:39` `CensusCubeQuery` | one class fronts all cube reads |
| Core method | `query.py:124-186` `_query()` | `cube.query(cond=, use_arrow=True, attrs=, dims=).df[dim_slice]` → pandas `DataFrame` |
| Criteria | `criteria.py` | `BaseQueryCriteria`, `CensusCubeQueryCriteria`, `MarkerGeneQueryCriteria` |
| Projection | `CensusCubeQueryParams` | column/dim projection — **reuse as-is** for column parity |
| Precedent backends | `query_duckdb.py`, `query_lance.py` | already implement this exact seam |

**8 public methods to reimplement per backend:**
`expression_summary`, `expression_summary_default`, `cell_counts`, `marker_genes`,
`cell_counts_df`, `cell_counts_diffexp_df`, `expression_summary_and_cell_counts_diffexp`
(returns a **tuple**, **bypasses `_query`**, slices by integer `group_id`), and the primary-dim
enumerators `list_primary_filter_dimension_term_ids` /
`list_grouped_primary_filter_dimensions_term_ids`.

**All consumers of the seam:**
- WMG API `backend/wmg/api/v2.py` — `query()` (L87/89/92), `markers()` (L183).
- DE API `backend/de/api/v1.py` — `differentialExpression()` (L253-263, L333), `filters()` via
  `cell_counts_df`.
- Pipeline `backend/wmg/pipeline/primary_filter_dimensions.py:52-61`.
- Note: WMG `filters()` reads snapshot dicts directly (not the seam) — unaffected.

The two return shapes any backend must satisfy: **(a)** categorical-dim result frames
(`expression_summary`/`cell_counts`/`marker_genes`), and **(b)** the integer-`group_id` diffexp
tuple. The API response builders downstream (`get_dot_plot_data` → `rollup`) are unchanged.

### Layer 2 — Write / build pipeline (`backend/wmg/pipeline/`)

Entrypoint/DAG: `pipeline.py:46-80` — 9 steps producing **7 TileDB arrays + 5 JSON sidecars**.

| # | Writer file | Array | Indexed dims | Key attrs | Schema |
|---|---|---|---|---|---|
| 1 | `expression_summary.py` | `expression_summary` | gene, tissue, organism | 7 categorical + nnz/sum/sqsum | `schemas/cube_schema.py:37` |
| 2 | `cell_counts.py` | `cell_counts` | tissue, organism | 7 categorical + n_cells | `cube_schema.py:71` |
| 3 | `expression_summary_default.py` | `expression_summary_default` | gene, tissue, organism | cell_type + nnz/sum/sqsum | `cube_schema_default.py` |
| 4 | `expression_summary_and_cell_counts_diffexp.py` | `expression_summary_diffexp` | `group_id` (uint32) | gene + sum/sqsum | `cube_schema_diffexp.py` |
| 5 | ″ (same file) | `expression_summary_diffexp_simple` | `group_id` | gene + sum/sqsum | `cube_schema_diffexp.py` |
| 6 | ″ (same file) | `cell_counts_diffexp` | cell_type, tissue, organism | n_cells + group_id/group_id_simple | `cube_schema_diffexp.py` |
| 7 | `marker_genes.py` | `marker_genes` | tissue, cell_type, organism | gene + marker_score/specificity | `marker_gene_cube_schema.py` |

**Load-bearing schema properties to reproduce (by layout, not just "an index"):** sparse arrays,
**row-major cell/tile order sorted on the indexed dims**, capacity=10000 tiling, dict+zstd-19 on
categoricals, byteshuffle+zstd-5 on numerics (`schemas/tiledb_filters.py`). Creation util:
`pipeline/utils.py:50` `create_empty_cube_if_needed`.

**Source read (out of scope — stays TileDB):** `expression_summary_and_cell_counts.py:50`
`open_soma(...)` + `axis_query("RNA", value_filter=...)`; obs/var pulled and renamed per
`constants.py`. The 5 JSON sidecars (`primary_filter_dimensions`, `dataset_metadata`,
`cell_type_ancestors`, `filter_relationships`, `cell_type_orderings`) are format-agnostic and
carry over unchanged.

### Layer 3 — Snapshot / load / deploy

| Concern | Location | Mechanism |
|---|---|---|
| Open/cache | `snapshot.py` `load_snapshot()` L115 → `_open_cube()` L394 | `tiledb.open`; global `cached_snapshot` L111; 7 arrays held open for process life; `cell_counts_df` + `cell_counts_diffexp_df` materialized to pandas at load (L342-343) |
| Versioning | `_get_latest_snapshot_id()` L185 | `latest_snapshot_identifier` text file at `s3://{WMG_BUCKET}/snapshots/v5/`; schema `v5` pinned in `constants.py:18`, `wmg/api/config.py:15`, `container_init.sh:15` |
| Build→S3 | `pipeline/load_cube.py` `upload_artifacts_to_s3()` | `aws s3 sync` **dir-to-S3** (retains 2 latest) |
| Runtime sync | `wmg/server/container_init.sh` | `aws s3 sync` → local cache `/single-cell-data-portal/census_cube_snapshot_cache` (~54 GB); else TileDB-VFS direct-S3 fallback |
| Runtime config | `data/tiledb.py` `create_ctx()` | `sm.tile_cache_size=50% VM`, `vfs.s3.*`, buffer sizes — **replacement rewrites this whole file** |

**Deploy constraints that bound the engine choice:** ECS Fargate, gunicorn + gevent; prod
16 GB / 4 vCPU / 5 workers / 3 instances; 100 GB ephemeral (54 GB cache); 60s healthcheck window
sized for the startup sync. An **embedded** engine must fit within worker memory and not block the
gevent event loop on I/O; a **server** engine breaks the dir-to-S3 snapshot model (needs an
always-on service + DB restore) — the biggest scope delta.

---

## 4. Options — weighing all candidates

The §2 requirement (layout must cluster the leading dims → result-size access) is the filter.

| Option | Status | Verdict |
|---|---|---|
| DuckDB / Parquet-as-table | **spiked** | 35–150× slower (full scan) — **OUT** |
| Lance (unsorted + sorted/clustered) | **spiked** | 1.2–10.7×, ~35 ms fixed floor, fails gate — **OUT** (closest) |
| chDB querying raw Parquet | analyzed | same scan trap (chDB only wins on its own MergeTree) — **OUT** |
| Zarr | analyzed | no query engine at all — **OUT** |
| In-memory pandas | analyzed | O(n) mask, scales with table size — **OUT** |
| **ClickHouse MergeTree — embedded via chDB (persistent MergeTree DB)** | **spiked 2026-07-02 — PASSED** | **the viable exit target** |
| ClickHouse MergeTree — server | not spiked | works, but breaks the dir-to-S3 model (always-on service); embedded is preferred |
| DIY: partition-by-gene + manifest index | analyzed | reimplements TileDB tiling in app code — only if dropping the C++ dep is the *sole* goal |

**Why chDB-embedded is the pick.** `ORDER BY (organism, gene, tissue)` builds a sparse primary
index; `WHERE gene IN (...)` binary-searches in-memory index marks and reads only the matching
8192-row granules — **result-size access**, the property §2 requires, and structurally closer to
TileDB's coordinate model than Lance's row-id/bitmap space. It keeps today's architecture:
in-process (no server), dir-to-S3 snapshot, embedded engine — the lazy way to get granule-skipping
without running a ClickHouse cluster. `sum`/`sqsum`/`nnz` are already precomputed, so a query is a
projected filter. It also aligns with the same SQL/ML-native ecosystem pull behind the Explorer Zarr work.

**Spiked on the real cube (2026-07-02) — PASSED the gate. First candidate that does.** Built
`ChdbCensusCubeQuery` behind the seam (`query_chdb.py`) + `wmg_build_chdb.py` (Parquet → MergeTree),
exported the real v5 staging snapshot (`1760292291`, `expression_summary` = 351,721,698 rows), and
ran the same real-cube harness Lance/DuckDB failed. **Parity: PASS on every shape** (main +
diffexp).

| shape | rows | TileDB | chDB | vs TileDB |
|---|---|---|---|---|
| default — genes only (60M) | 1,430 | 7.9 ms | 4.9 ms | **1.6× faster** |
| genes + tissues (351M) | 90 | 13.9 ms | 3.9 ms | **3.6× faster** |
| genes + cell_types (secondary, 351M) | 103 | 8.5 ms | 10.0 ms | 1.2× slower |
| diffexp — full group_ids | 39,372 | 19.6 ms | 10.2 ms | **1.9× faster** |
| diffexp — simple group_ids | 39,372 | 16.1 ms | 6.4 ms | **2.5× faster** |

Granule diagnostic (the exact test that killed Lance): gene-selective lookups prune to **3 / 42,969
granules** in ~2–3 ms. Where Lance's scalar index returned row-ids but then did a scattered take with
a ~35 ms floor, chDB's MergeTree primary key does a contiguous granule slice — the same result-size
access TileDB gets from coordinate tiling. This is the layout property §2 demanded, finally
reproduced.

**Caveats (honest bounds on the result).** (1) The secondary `cell_type` filter is the one shape
~1.2× slower — marginal and absolute-fine (10 ms), tunable via a `set` index or adding cell_type to
the sort key. (2) DB is ~17 GB vs TileDB ~7.6 GB (~2.2×). (3) Measured warm-cache, single-process,
local disk — same conditions as the TileDB baseline (fair), but **not yet tested under gevent
concurrency or at deploy scale**; that's the next thing to verify before committing to the exit.
(4) The `organism`-only diagnostic number is a `count()` optimization, not a data-materializing
measure — the main-table cases are the honest latency.

---

## 5. Scope of rework

**Query (~1 file).** Reimplement the 8 methods behind the seam for one backend, mirroring
`query_lance.py`. Must cover both shapes: categorical-dim `expression_summary` **and** the integer
`group_id` diffexp path. Reuse `CensusCubeQueryParams` unchanged for column parity.

**Write (7 writers).** Replace the 7 array writers with new-format writers. The critical
constraint: **preserve the sort/cluster layout on the indexed dims** — that's the load-bearing
property, not incidental. Keep the 5 JSON sidecars and the source read as-is.

**Snapshot/deploy.** Rewrite `_open_cube` and `create_ctx` (`data/tiledb.py`) for the new engine;
keep the `load_snapshot` caching shape and the versioning contract. With chDB-embedded, keep the
`aws s3 sync` dir model — **no terraform change**. A server engine would need new infra
(service, DB restore, healthchecks) — treat as a separate, larger track.

**Unchanged:** criteria classes, `CensusCubeQueryParams`, the 5 JSON sidecars, the source census
read, the API response builders (`get_dot_plot_data`/`rollup`), and the snapshot versioning
contract (`latest_snapshot_identifier`, schema-version pinning).

---

## 6. Benchmarks & the hard gate

**Done — the harness ran and chDB passed (2026-07-02).** Reused the `scripts/wmg_*_spike_realcube.py`
pattern on the **real v5 staging snapshot** (`1760292291`, `expression_summary` = 351,721,698 rows):
`query_chdb.py` behind the seam, `wmg_build_chdb.py` (TileDB → Parquet → MergeTree), the realcube
runner, and `wmg_chdb_spike_diffexp.py` for the diffexp path. Full results table in §4.

**The gate (what flipped DuckDB from "faster" to "150× slower"):** median latency for all query
shapes must be **≤ TileDB's**. Result: **chDB met it** — faster on 4 of 5 shapes, 1.2× slower on the
secondary `cell_type` shape only.

| Shape | Cube size | TileDB baseline | chDB | gate |
|---|---|---|---|---|
| default — genes only | 60M rows | 6–10 ms | 4.9 ms | ✅ |
| genes + tissues | 351M rows | 9–14 ms | 3.9 ms | ✅ |
| genes + cell_types | 351M rows | 6–15 ms | 10.0 ms | ⚠️ 1.2× (marginal) |
| diffexp (full / simple) | 254M / 61M rows | 16–20 ms | 6–10 ms | ✅ |

The **per-predicate diagnostic** that exposed Lance's gap confirmed the primary-key prefix prunes to
3 / 42,969 granules for gene-selective lookups (§4). Concurrency + deploy-scale were also run — see §9;
they change *how you deploy* chDB, not whether the latency holds.

---

## 7. Testing & parity

**Parity:** chDB is numerically identical to TileDB for `expression_summary`
(default/indexed/secondary-filter variants), `cell_counts`, and the diffexp tuple — **verified PASS
on the real cube** (the diffexp check needed the spurious `reset_index()` "index" column dropped; the
underlying data was already identical). `marker_genes` parity not separately asserted in the run but
uses the same `_query` path. This is the same parity bar DuckDB and Lance also passed — the
difference is chDB also passes on **latency**, which they didn't.

**Test surface:** the seam contract itself — the 8 methods and their DataFrame shapes. Existing
unit tests under WMG/DE that construct `CensusCubeQuery` exercise this; a new backend must pass
them unchanged. The seam being one class is what makes the swap testable in isolation (the
DuckDB/Lance spikes are the proof).

---

## 8. Residual: the source read stays TileDB

Even a perfect cube exit does **not** remove TileDB from the deployment, because the pipeline
ingests the Census corpus via `tiledbsoma`. Census publishes one TileDB-free artifact — the
per-dataset source H5ADs — but those are *pre-integration* raw inputs; reproducing the pipeline's
cross-corpus `axis_query` over them means re-harmonizing every dataset yourself (a fork of CZI's
integration, cf. `chanzuckerberg/multimodal-slicing`'s gene-universe work noted in the findings
doc). The read leaves TileDB only if (a) CZI republishes the integrated corpus in a non-TileDB
format upstream, or (b) that fork is undertaken. Both are separate, larger efforts outside this
cube-exit scope.

---

## 9. Concurrency & deployment architecture

The realcube latency numbers (§4) were single-query, single-process. Prod is 5 gunicorn/gevent workers
in a 16 GB Fargate task. A concurrency spike (`scripts/wmg_chdb_spike_concurrency.py`, process pool +
per-worker DB, RSS sampling, gevent check) surfaced a hard architectural fact:
**embedded chDB is the wrong *serving* engine under concurrency — but the MergeTree it builds is not.**

**What the spike measured (real cube, this laptop):**

| backend | workers | thru q/s | p50 | p95 | p99 | peak RSS |
|---|---|---|---|---|---|---|
| TileDB | 1 / 5 / 10 | 83 / 230 / 185 | 6 / 9 / 17 ms | 11 / 18 / 37 ms | 18 / 26 / 58 ms | 0.2 / 1.4 / 2.5 GB |
| chDB | 1 | 175 | 4.1 ms | 5.6 ms | 7.6 ms | 0.8 GB |
| chDB | 3+ | — | — | — | — | **fails** (see below) |

Three findings, none about latency:
1. **chDB takes an exclusive lock on its data dir** (`Cannot lock file …/status`). Multiple processes
   **cannot share one on-disk DB** — the "5 workers share the S3-synced snapshot" model doesn't work.
2. **Per-worker copies are wasteful *and* fragile.** N workers → N × ~17 GB disk, and cloning a
   persisted chDB dir fails to reopen (`Directory metadata/wmg already exists` / `EmbeddedServer
   BAD_ARGUMENTS`) — you'd have to rebuild each copy from Parquet.
3. **In-process concurrency serializes.** 10 gevent greenlets took **14.7×** a single warm query — a
   native chDB call blocks the event loop, so a worker gets zero query overlap. (TileDB shares a
   read-only dir fine and scales to 5 workers; that's the incumbent model chDB-embedded can't match.)

**The fix (ClickHouse-native, confirmed in CH docs): build embedded, serve with a server sidecar over
a read-only disk.**

```
Pipeline (offline):  build MergeTree with chDB/clickhouse-local → clickhouse-static-files-uploader
                     → s3://…/snapshots/vN/<id>/   (mirrors today's dir-to-S3 + latest_snapshot_identifier)
ECS task:
  ├─ gunicorn+gevent workers — CensusCubeQuery talks to CH over localhost (HTTP/native client)
  └─ clickhouse-server sidecar — read-only `web` / `s3_plain` disk over the S3 snapshot + local cache disk
```

ClickHouse's **read-only `web`/`s3_plain` disk** is purpose-built for this: *"a read-only disk, its
data is only read and never modified,"* prepared offline with `clickhouse-static-files-uploader`, and
**multiple independent instances can attach the same data concurrently** (immutable → no consistency
issue), with an LRU **local cache disk** for hot granules. This resolves every finding at once:

| Spike finding | Sidecar + read-only disk |
|---|---|
| dir lock (can't share a dir) | one server owns the data; workers are clients |
| per-worker memory ×N | one server, **shared** caches — not multiplied |
| per-worker disk copies (N×17 GB) | read-only S3 disk + local cache — no copies |
| **gevent serialization (14.7×)** | queries become **localhost socket calls** → gevent yields; the in-process blocking disappears |
| OOM crash under load | server enforces `max_concurrent_queries` + memory limits → backpressure |

**Cost:** it's a localhost **sidecar container** in the ECS task def — not pure-embedded, but not a
managed cluster either (no separate service, no replication). Two variants: (a) `web` disk straight
from S3 + local cache (fast cold start, S3 latency on cold granules); (b) sync the static MergeTree to
local disk as `container_init.sh` does today, serve read-only over the local path (local-disk latency +
server concurrency — closest to current ops). Default to (b). chDB stays the **build-side** engine.

---

## 10. Recommendation

**chDB-embedded answers the hard question; the deployment shape is a clickhouse-server sidecar.** The
read-path risk this whole effort circled — "can any non-TileDB engine meet the selective-lookup
latency?" — is answered: **yes**, parity on every shape and faster than TileDB on 4 of 5, because the
MergeTree primary index reproduces TileDB's result-size granule access. Concurrency then showed the
*serving* layer must be a server sidecar over a read-only disk (§9), not the embedded engine in each
worker. Both risky questions are now de-risked.

**What's left is real work, not research** (in order):
1. **Full rework** per §5 — the 7 writers, snapshot wiring, and swap `_open_cube`/`create_ctx` for a CH
   client. The query seam (`query_chdb.py`) is already written and validated.
2. **Stand up the sidecar** — clickhouse-server in the ECS task with a read-only `web`/`s3_plain` (or
   local read-only) disk; re-run the concurrency harness against it (queries as socket calls) to
   confirm p99 under load and shared-cache memory fit the 16 GB task.
3. **Tune the secondary `cell_type` shape** if the 1.2× matters (`set` skip-index or sort-key change),
   and `LowCardinality(String)` on the ontology-ID columns to shrink the ~17 GB DB.

**But still: don't start that work absent a driver.** The cube meets its budget on TileDB today, and a
cube exit doesn't remove TileDB from the deployment anyway (the source read is upstream SOMA, §8). The
value of this spike is *optionality* — if a driver lands (dropping the C++ dep, an object-store
mandate), the storage-and-serving path is now known end-to-end, not a research project. Keep TileDB
now; reach for the chDB-build + CH-sidecar path when forced.

---

_Companion spike artifacts (throwaway, not production): `query_chdb.py` (the one that passed),
`query_duckdb.py`, `query_lance.py`, `scripts/wmg_build_chdb.py`, `scripts/wmg_chdb_spike_realcube.py`,
`scripts/wmg_chdb_spike_diffexp.py`, `scripts/wmg_chdb_spike_concurrency.py` (process-pool load + RSS +
gevent check), and the Lance/DuckDB equivalents. See the findings doc's Artifacts section._
