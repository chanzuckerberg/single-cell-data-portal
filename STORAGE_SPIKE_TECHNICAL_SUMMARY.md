# Storage-backend spike — technical summary

The single technical synthesis of the TileDB spike: the question, the method, the results, and the
path forward per component. Deeper detail lives in three focused docs:
[WMG cube exit plan](WMG_CUBE_TILEDB_EXIT_REWORK.md) ·
[engine architecture](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md) ·
[census corpus scope](CENSUS_CORPUS_GENERATION_SCOPE.md). For the non-technical version see the
[product brief](STORAGE_SPIKE_PRODUCT_BRIEF.md).

Branch: `spike/tiledb-storage-backend-migration` (pushed to `chanzuckerberg/single-cell-data-portal`,
no PR, not merged). All three repos in scope (`single-cell-data-portal`, `cellxgene-census`,
`single-cell-explorer`) are CZI — boundaries here are component/repo/team, not organizational.

**Question:** can the single-cell stack move off TileDB, and how? **Answer:** it depends on the data's
*geometry*, and each component has a proven or validated path. This is a **de-risking spike** that
maps the options and de-risks the path — whether and when to move is a product/leadership decision this
analysis informs, not one it makes.

---

## The three artifacts — match the format to the geometry

The stack isn't one storage problem. Three artifacts have different shapes, so each wants a different
backend, even though all three use TileDB today.

| Artifact | Repo (CZI) | Geometry | Finding | Deep dive |
|---|---|---|---|---|
| **Explorer per-dataset matrix (`.cxg`)** | single-cell-explorer | per-cell **tensor** | **Zarr works** — proven, ran in rdev | below |
| **WMG cube** | single-cell-data-portal | sparse predicate-filtered **OLAP aggregate** | **chDB validated (matches TileDB); DuckDB/Lance/Zarr regress; not yet built** | [cube plan](WMG_CUBE_TILEDB_EXIT_REWORK.md) |
| **Census corpus** (the cube's source) | cellxgene-census | integrated per-cell **tensor** | **Feasible; heaviest** — datacenter-scale build + public API | [census scope](CENSUS_CORPUS_GENERATION_SCOPE.md) |

Two independent axes underlie this: **tensor vs OLAP-aggregate** (the geometry) and **TileDB-SOMA vs
plain TileDB** (the API). The `.cxg` is tensor-shaped on *plain* TileDB; the corpus is tensor-shaped on
*SOMA*; the cube is OLAP on plain TileDB. Tensors (`.cxg`, corpus) fit Zarr; the OLAP cube does not.

---

## The gate — parity + latency on the real cube

The decisive test: run every query shape against the **real v5 staging snapshot** (`1760292291`,
`expression_summary` = 351,721,698 rows) and require median latency **≤ TileDB's** (6–15 ms). This is
what flipped DuckDB from "2× faster" on a toy fixture to "150× slower" on the real cube — nothing ships
on theory. The cube's runtime pattern is selective point-lookups returning ~1–6k rows out of 352M,
under a <10 s interactive budget.

## Results

| Candidate | Result vs TileDB | Root cause |
|---|---|---|
| **Zarr** (for the cube) | n/a | no query engine — fine for tensors, wrong for the OLAP cube |
| **DuckDB / Parquet** | ❌ 35–150× slower | scan-first; ART index unused for `IN`-lists |
| **Lance** (sorted + unsorted) | ❌ 1.2–10.7× slower, ~35 ms floor | scalar index returns row-ids, data not clustered → scattered take |
| **chDB** (embedded ClickHouse MergeTree) | ✅ **parity everywhere; faster on 4 of 5 shapes** | sorted MergeTree + sparse primary index → contiguous granule slice |

chDB is the first backend whose storage layout reproduces what makes TileDB fast: a gene-selective
lookup **prunes to 3 of 42,969 granules** — result-size access, not a scan. Storage is ~17 GB vs
TileDB's ~7.6 GB (tunable). Full per-shape latencies and the granule diagnostic are in the
[cube exit plan](WMG_CUBE_TILEDB_EXIT_REWORK.md).

**Why chDB and not "a columnar DB":** the deciding property is the storage **layout** (does it cluster
the leading dims so a lookup is a contiguous slice?), not the "columnar" label — which is exactly why
two *columnar* candidates (DuckDB, Lance) lost. See [engine architecture](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md).

## Concurrency → a server sidecar

The latency numbers were single-process. Prod is 5 gunicorn/gevent workers in a 16 GB task, and
**embedded chDB is the wrong serving engine** there (exclusive data-dir lock; 10 greenlets serialize
14.7×). The fix, validated locally: build with embedded chDB → publish a read-only MergeTree to S3 →
serve with a **clickhouse-server sidecar** over a read-only `web`/`s3_plain` disk. Then gevent overlaps
(4.1×, not 14.7×), one ~296 MB server serves all workers, and throughput (280–327 q/s) beats TileDB's
5-worker 230 q/s. Detail in the [cube exit plan](WMG_CUBE_TILEDB_EXIT_REWORK.md).

## Explorer `.cxg` → Zarr (proven)

The Explorer per-dataset matrix is a dense tensor with positional slicing and no query engine — a clean
Zarr fit. Proven in the separate `single-cell-explorer` repo (branch `poc/zarr-adaptor`, PR #1369, ran
in rdev). CZI's `chanzuckerberg/multimodal-slicing` spike independently corroborates the same
tensor→Zarr thesis on census-derived per-dataset H5ADs.

## The residual — the Census corpus

Even a perfect cube exit leaves TileDB in the offline build, because the cube's source is the Census
corpus, read via `open_soma` → `tiledbsoma`. That corpus is a 512-GiB weekly build in the
`cellxgene-census` repo, and the same TileDB-SOMA object is the storage engine under the **public
`cellxgene-census` reader API** — read live over S3 by a broad external community. So its format is a
stable public contract and the heaviest, most coordination-sensitive component to move. Full treatment,
including the only external-facing TileDB dependency, in the [census corpus scope](CENSUS_CORPUS_GENERATION_SCOPE.md).

---

## Path forward (per component)

Each component has a proven or validated path off TileDB — a move is execution, not research. Whether
and when to move is a product/leadership decision.

- **Explorer `.cxg`** → **Zarr** (proven).
- **WMG cube** → build with embedded chDB → publish a read-only MergeTree to S3 → serve with a
  clickhouse-server sidecar. What's left is implementation (7 pipeline writers, snapshot wiring, swap
  the `CensusCubeQuery` seam) plus two things a laptop can't verify: real-S3 cold-read latency and
  behavior on ECS hardware.
- **Census corpus** → the heaviest, most coordination-sensitive component. Paths (lowest-effort first):
  a non-TileDB SOMA backend in the Census stack, or a source-H5AD bypass — not a data-portal fork of
  the `cellxgene-census` builder. One CZI program across the three components.

Note: removing TileDB *entirely* requires both the cube exit *and* the corpus change — either alone
leaves the C++ core via the other.

---

## Artifacts & repro

**Backend seam (throwaway spikes):** `backend/common/census_cube/data/query_{chdb,duckdb,lance}.py` —
all behind the single `CensusCubeQuery._query` seam that bounds the swap.

**Scripts (`scripts/`):** `wmg_build_chdb.py`, `wmg_chdb_spike_realcube.py`, `wmg_chdb_spike_diffexp.py`,
`wmg_chdb_spike_concurrency.py` (embedded), `wmg_ch_sidecar_concurrency.py` (server sidecar); plus the
DuckDB/Lance equivalents and `wmg_export_cube_to_parquet.py`.

**Repro needs:** the spike venv (`venv-311/bin/python`), a local snapshot (interactive Okta login for
S3), `pip install chdb clickhouse-connect`, single-binary ClickHouse (`curl https://clickhouse.com/ | sh`).
Explorer `.cxg`→Zarr lives in `single-cell-explorer` (branch `poc/zarr-adaptor`, PR #1369).
