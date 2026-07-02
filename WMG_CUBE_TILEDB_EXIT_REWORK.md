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
| **ClickHouse MergeTree — embedded via chDB (persistent MergeTree DB)** | **not spiked** | **RECOMMENDED path to spike if forced** |
| ClickHouse MergeTree — server | not spiked | works, but breaks the dir-to-S3 model (always-on service) |
| DIY: partition-by-gene + manifest index | analyzed | reimplements TileDB tiling in app code — only if dropping the C++ dep is the *sole* goal |

**Why chDB-embedded is the pick.** `ORDER BY (organism, tissue, gene)` builds a sparse primary
index; `WHERE gene IN (...)` binary-searches in-memory index marks and reads only the matching
8192-row granules — **result-size access**, the property §2 requires, and structurally closer to
TileDB's coordinate model than Lance's row-id/bitmap space. It keeps today's architecture:
in-process (no server), dir-to-S3 snapshot, embedded engine — the lazy way to get granule-skipping
without running a ClickHouse cluster. `sum`/`sqsum`/`nnz` are already precomputed, so a query is a
projected filter, and server-side aggregation could shrink payloads. It also aligns with the same
SQL/ML-native ecosystem pull behind the Explorer Zarr work.

**Why it's still unproven.** This is exactly the "looks right in theory" position DuckDB was in
before benchmarking flipped it from "2× faster" (toy fixture) to "150× slower" (real cube). The
load-bearing risk to verify: that an `IN`-list on the **primary-key prefix** actually prunes
granules on the real 351M-row snapshot (not a full-part scan), and that the secondary
low-cardinality filters don't reintroduce Lance's bitmap-intersection cost. Nothing ships on
theory (§6).

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

**Reuse the spike harness** — the `scripts/wmg_*_spike_realcube.py` pattern that already exists for
DuckDB and Lance: parity + latency on the **real v5 staging snapshot** (`1760292291`,
`expression_summary` = 351,721,698 rows). Adding a chDB backend is a new `query_chdb.py` behind the
seam plus a `wmg_build_chdb.py` exporter (TileDB → Parquet → MergeTree), then the same runner.

**The gate (what flipped DuckDB from "faster" to "150× slower"):** median latency for all three
query shapes must be **≤ TileDB's**:

| Shape | Cube size | TileDB baseline |
|---|---|---|
| default — genes only | 60M rows | 6–10 ms |
| genes + tissues | 351M rows | 9–14 ms |
| genes + cell_types | 351M rows | 6–15 ms |

Plus the **per-predicate diagnostic** that exposed Lance's structural gap: measure `gene`-only,
`organism`-only, and the AND-combinations separately, to confirm the primary-key prefix prunes and
low-cardinality dims don't blow up. If any shape regresses, the candidate is out — same bar Lance
failed.

---

## 7. Testing & parity

**Parity:** numerically identical to TileDB for `expression_summary`
(default/indexed/secondary-filter variants), `cell_counts`, `marker_genes`, and the diffexp tuple
— the same parity check the DuckDB and Lance spikes passed (PASS everywhere; latency was the only
failure).

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

## 9. Recommendation

**Don't spike now.** No driver forces a TileDB exit, the cube meets its latency budget on TileDB,
and the source read stays TileDB regardless — so a cube exit doesn't even remove the dependency
from the deployment. This doc is the runbook, not a work order.

**If an exit is forced:** benchmark **chDB-embedded first** — it's the only untested candidate
whose storage model can plausibly meet the ≤-TileDB gate, and it preserves the embedded,
dir-to-S3, no-server architecture. Run it through the existing real-cube harness and the hard gate
in §6 before writing any writer or deploy code. If chDB fails the gate the way Lance did, the
honest fallback is DIY partition-by-gene — i.e. rebuilding TileDB's tiling yourself — at which
point keeping TileDB is the lazy and correct answer.

---

_Companion spike artifacts (throwaway, not production): `query_duckdb.py`, `query_lance.py`,
`scripts/wmg_*_spike_realcube.py`, `scripts/wmg_build_lance*.py`. See the findings doc's Artifacts
section._
