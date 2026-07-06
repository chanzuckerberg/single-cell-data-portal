# Storage-backend spike — overview (why / what / recommendation)

The one-page entry point to the TileDB storage spike. Distills the whole effort; each claim links
to the detailed doc that backs it. Read this first; go deeper only where you need to.

**Detailed docs:** [investigation summary](STORAGE_SPIKE_SUMMARY.md) ·
[findings & evidence](STORAGE_BACKEND_MIGRATION_FINDINGS.md) ·
[cube-exit rework](WMG_CUBE_TILEDB_EXIT_REWORK.md) ·
[ClickHouse vs TileDB architecture](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md) ·
[census corpus generation](CENSUS_CORPUS_GENERATION_SCOPE.md)

---

## Why

**Question: can the single-cell stack move off TileDB?** TileDB is a C++ storage engine used across
three CZI components — the WMG cube and Explorer matrices (in `single-cell-data-portal` and
`single-cell-explorer`) and the Census corpus they read (in `cellxgene-census`). *(All three are CZI
`chanzuckerberg` repos; the boundaries here are component/team, not organizational.)* A move would
only be driven by something specific — dropping the native C++ dependency, or an org-wide
object-store / ML-native (Zarr/Arrow/SQL) storage mandate. **Nothing forces that today**, so this was
a **de-risking spike, not a migration**: prove whether an exit is possible and know the path, before
any driver lands.

## What we found

The stack isn't one storage problem — it's **three artifacts with different geometry**, so each gets
its own answer. Full evidence in the [findings doc](STORAGE_BACKEND_MIGRATION_FINDINGS.md).

| Artifact | Geometry | Verdict | Detail |
|---|---|---|---|
| **Explorer per-dataset matrix (`.cxg`)** | per-cell **tensor** | **Zarr works** — proven in `single-cell-explorer` (PR #1369, ran in rdev) | [findings §1](STORAGE_BACKEND_MIGRATION_FINDINGS.md) |
| **WMG cube** | sparse predicate-filtered **OLAP aggregate** | **Keep TileDB now; a chDB exit is viable if forced** | [rework](WMG_CUBE_TILEDB_EXIT_REWORK.md) · [architecture](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md) |
| **Census corpus (source)** | integrated per-cell **tensor** | **Leave as-is** — separate CZI repo (`cellxgene-census`) + a public API; the data-portal only reads it offline | [census scope](CENSUS_CORPUS_GENERATION_SCOPE.md) |

The load-bearing lesson from the cube work: **the property that matters is the storage *layout*, not
the query engine or the "columnar" label.** TileDB is fast because it clusters/tiles data on the
indexed dims, so a selective lookup is a contiguous slice (result-size access), not a scan. Any
replacement must reproduce that. Benchmarked on the **real 351M-row cube** against a hard gate
(median latency ≤ TileDB's, 6–15 ms):

- DuckDB/Parquet — **35–150× slower** (scan-first). ❌
- Lance — **1.2–10.7× slower**, ~35 ms floor (row-id/bitmap, not clustered). ❌
- **chDB (embedded ClickHouse MergeTree) — PASSED**, faster than TileDB on 4 of 5 shapes ✅
  (its sparse primary index prunes to 3 of 42,969 granules — the same result-size access).

Concurrency then reshaped the *deployment*, not the verdict: embedded chDB can't serve 5 gevent
workers (dir lock, event-loop blocking), but a **clickhouse-server sidecar over a read-only S3 disk**
fixes it (validated locally). See the [architecture doc](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md) for
why the two engines behave differently, and the [rework doc](WMG_CUBE_TILEDB_EXIT_REWORK.md) for the
layer-by-layer exit scope.

**The residual:** even a perfect cube exit leaves TileDB in the offline build, because the cube's
source is the Census corpus (TileDB-SOMA, read via `open_soma`). That corpus is a 512-GiB weekly build
in a different CZI repo (`cellxgene-census`) — and, though the data-portal's read of it is offline
(pipeline-only, never on the data-portal serving path), the same TileDB-SOMA object is also the storage
engine under the **public `cellxgene-census` reader API** (Python/R, read live over S3 by a broad
community that includes external, non-CZI users). So the corpus format is a stable public contract
owned by another CZI team — reinforcing "leave it," and making a non-TileDB **SOMA backend** (in the
Census stack) the only clean exit. Scoped in
[census corpus generation](CENSUS_CORPUS_GENERATION_SCOPE.md).

## Recommendation

**Keep TileDB today.** No driver forces a change, every artifact meets its budget, and the value of
the spike is **optionality** — the exit path per artifact is now known and de-risked, not a research
project:

- **Explorer `.cxg`** → **Zarr** (proven).
- **WMG cube** (only if forced) → **build with embedded chDB → publish a read-only MergeTree to S3 →
  serve with a clickhouse-server sidecar** over a read-only `web`/`s3_plain` disk. What's left is
  implementation (7 pipeline writers, snapshot wiring, swap the `CensusCubeQuery` seam) plus two
  things a laptop can't verify: real-S3 cold-read latency and behavior on ECS hardware.
- **Census corpus** → **leave as TileDB-SOMA.** If a full TileDB purge is ever mandated, the
  lowest-effort paths are a non-TileDB SOMA backend in the Census stack or a source-H5AD bypass —
  **not** a data-portal fork of the `cellxgene-census` builder. This is one CZI program across the
  three components, led with the `cellxgene-census` team, not an external ask.

**In one line:** don't migrate now; reach for these paths only when a concrete driver appears.

---

_Spike branch: `spike/wmg-parquet-duckdb` (pushed to `chanzuckerberg/single-cell-data-portal`, no PR,
not merged). Throwaway artifacts (spike query backends, build/benchmark scripts) are inventoried in
the [investigation summary](STORAGE_SPIKE_SUMMARY.md) and [findings doc](STORAGE_BACKEND_MIGRATION_FINDINGS.md)._
