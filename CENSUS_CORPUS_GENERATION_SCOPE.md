# Census corpus generation — scoping the storage-spike residual

Sibling to [`STORAGE_SPIKE_SUMMARY.md`](STORAGE_SPIKE_SUMMARY.md) (the spike verdict),
[`WMG_CUBE_TILEDB_EXIT_REWORK.md`](WMG_CUBE_TILEDB_EXIT_REWORK.md) (the cube-exit rework), and
[`CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md`](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md) (the engine
comparison). Every one of those docs stops at the same unresolved edge — the WMG cube's *source*
is the Census corpus, which is TileDB. This doc scopes that edge: **how the Census corpus is
generated, where its TileDB dependency actually sits, and what moving it off TileDB would take.**

**This is a decision doc, not a spike.** A real Census build needs **512 GiB RAM / 2 TiB swap /
1.8 TiB disk** (§2), so unlike the cube work there is no laptop-runnable proof here — the analysis
is on the code and the published build requirements, not benchmarks.

---

## 1. Context & scope

The cube spike proved the WMG **cube** can leave TileDB (chDB build + clickhouse-server sidecar).
But both `WMG_CUBE_TILEDB_EXIT_REWORK.md` §8 and `STORAGE_SPIKE_SUMMARY.md` §7 close on the same
residual:

> the cube's *source* is the Census corpus, read via `cellxgene_census.open_soma` → `tiledbsoma`
> (TileDB). Even a perfect cube exit leaves TileDB in the pipeline.

**Corpus vs cube — the distinction this whole doc rests on.**

- The **cube** is WMG's own artifact: 7 precomputed TileDB arrays of *aggregate* stats
  (`nnz`/`sum`/`sqsum`, cell counts) keyed by gene/tissue/organism/cell-type/etc. Built and owned
  in *this* repo (`backend/wmg/pipeline/`). The cube spike is entirely about this.
- The **corpus** is the upstream CZI **Census**: the harmonized, integrated per-cell
  expression matrix for all of CELLxGENE, published as a TileDB-SOMA object. Built and owned in a
  **different, CZI-owned repo** — `../cellxgene-census/tools/cellxgene_census_builder/`. This repo
  only *reads* it.

The portal builds **no corpus of its own.** WMG opens the upstream Census read-only and aggregates
it down to the cube:

```
backend/wmg/pipeline/expression_summary_and_cell_counts.py:50
    cellxgene_census.open_soma(census_version="latest")   # CensusParameters.census_version
      → axis_query("RNA", obs_query=AxisQuery(value_filter=...))   # per organism
      → aggregate to nnz/sum/sqsum → write cube
```

The `value_filter` (`backend/wmg/pipeline/constants.py:130` `CensusParameters`) selects
`is_primary_data == True and nnz >= <threshold> and cell_type != 'unknown'`, minus a handful of
system-level tissues. (Note: the pipeline's `corpus_path` variable is a **misnomer** — it's the
local *cube output* directory, not an input corpus.)

So "moving the Census corpus off TileDB" is **not a change in this repo** — it's a change to an
external, CZI-owned build that runs weekly at datacenter scale. That framing bounds every option in §4.

---

## 2. How the Census corpus is generated today

All paths below are in the sibling repo `../cellxgene-census/`.

**Builder location & entry points** — `tools/cellxgene_census_builder/`:
- Orchestrator: `src/cellxgene_census_builder/__main__.py` (command chosen by the
  `CENSUS_BUILD_COMMAND` env var: `build`, `full-build`, `release`, …). `full-build` chains
  prebuild-checks → **build-soma** → reports → S3 copy → release → cleanup.
- The SOMA builder itself: `build_soma/build_soma.py::build()` — the 5-step driver, run inside a
  dask distributed client.
- Container: `entrypoint.sh` (pins `OMP/BLAS/MKL_NUM_THREADS=1`) → `python -m cellxgene_census_builder .`.

**Output format — TileDB-SOMA.** `tiledbsoma==1.11.4`, pinned in `pyproject.toml` and deliberately
lagged behind the reader `cellxgene-census` package so readers can always open builder output. The
object layout and all TileDB platform config (compression filters, tile sizes, capacity `2^16`,
row-major) live in `build_soma/globals.py`:

- Root `soma.Collection` with two sub-collections:
  - **`census_info`** — `datasets`, `summary`, `summary_cell_counts`, `organisms` dataframes.
  - **`census_data`** — one `soma.Experiment` per organism (`homo_sapiens`, `mus_musculus`), each with:
    - `obs` dataframe (cell metadata),
    - `ms/RNA` measurement holding `var` (gene metadata), the `X` layers **`raw`** and
      **`normalized`** (float32 sparse `SparseNDArray`), and a boolean
      `feature_dataset_presence_matrix`.
- Schema constants: `CENSUS_SCHEMA_VERSION = "2.1.0"`, ingesting CELLxGENE dataset schema
  `CXG_SCHEMA_VERSION = "5.2.0"`.

**The 5-step build** (`build_soma.py::build()`):
1. Get source datasets — load manifest from the CxG REST API
   (`.../curation/v1/datasets?schema_version=5.2.0`) or a CSV; apply `dataset_blocklist.txt`;
   download/stage all H5ADs in parallel (dask/fsspec).
2. Create root collection + empty experiments/measurements.
3. Populate obs & var axes — concat obs / union var across all H5ADs, assign `soma_joinid`s.
4. Populate X layers — write `raw` + `normalized` sparse matrices (parallel over datasets), then the
   presence matrix. **The memory-bound step.**
5. Save axis & summary info (obs/var dataframes, dataset manifest, summary cell counts).
   Then optional TileDB **consolidate + vacuum** and a full re-read **validate**.

**Dependencies:** `tiledbsoma`, `anndata` (H5AD read), `dask`+`distributed` (parallel build),
`pyarrow`/`pandas`/`numpy`/`scipy`/`numba`, `cellxgene-ontology-guide`, `s3fs`/`fsspec`. (`tiledb`
comes transitively via `tiledbsoma`.)

**Build scale** (`build_state.py` host-validation minimums + `README.md`): **512 GiB** physical
memory, **2 TiB** swap, **1.8 TiB** free disk; dask over up to ~48 workers; step 4 budgets **64 GiB
per worker** (`n_workers = total_RAM // 64 GiB`). Corpus size: ~74M human + ~41M mouse cells across
~800 datasets (LTS 2024-07-01).

**Release model** (`release_manifest.py`, `dev-docs/cellxgene_census_storage_and_release_policy.md`):
published to `s3://cellxgene-data-public/cell-census/<tag>/{soma,h5ads}/`; a `release.json` maps
aliases → releases (`latest` = weekly, ~1-month retention; `stable`/`V#` = LTS, 6-month cadence,
5-year support). WMG pins `census_version = "latest"`.

---

## 3. Where TileDB actually sits — build-time vs serve-time

The load-bearing distinction that shrinks this residual's urgency:

- **Serving** (the WMG/DE API tier) does **not** touch the corpus. It reads only the cube. So a
  TileDB-free *serving* tier is already achievable via the chDB cube — no corpus change needed.
- **The corpus read is offline, pipeline-only** — it happens once per snapshot build in the WMG
  pipeline job (`open_soma` above), not on any request path.

So the residual is a **build-image dependency**, not a serving one: `tiledbsoma` must be installed
in the offline cube-build job because that's what reads the Census. Removing it means the *upstream
corpus* must exist in a non-TileDB form.

**Why `open_soma` == TileDB, structurally.** SOMA is an **API spec** (stack-agnostic by design);
`tiledbsoma` is currently the **only production implementation**. So today "read the Census via SOMA"
and "read TileDB" are the same statement — not by WMG's choice but because no other SOMA backend
ships. That is exactly the seam option (d) in §4 leans on.

---

## 4. Options

The §1 framing — external, CZI-owned, 512-GiB weekly build — is the filter. "What it takes" is
measured in *who owns the work* and *how much*, not in latency (there's nothing to benchmark).

| # | Option | What it takes | Cost | Verdict |
|---|---|---|---|---|
| a | **Leave it — consume upstream TileDB-SOMA as-is** | Nothing. `open_soma` stays; TileDB lives only in the offline build image | ~0 | **Recommended** |
| b | **Bypass the corpus — build the cube from source H5ADs** | Ingest the pre-integration H5ADs Census already publishes (`.../<tag>/h5ads/`), then re-implement CZI's cross-dataset **integration + gene-universe harmonization** to reproduce today's `axis_query` semantics | Large, **portal-owned, ongoing** (tracks upstream schema/integration changes forever) | Fallback only if forced |
| c | **Fork the builder to emit Zarr/Parquet** (obs/var → Parquet, X → sparse Zarr) | Fork `cellxgene_census_builder`, replace the SOMA/TileDB writer layer in `build_soma/`; still run the 512-GiB weekly build; WMG re-aggregates from the new format | **Largest.** Forking a datacenter-scale upstream build CZI owns and reruns weekly | Only as an *upstream contribution*, never a portal fork |
| d | **Wait for / contribute a non-TileDB SOMA backend** | If an Arrow/Zarr-backed SOMA implementation lands upstream, `open_soma` becomes TileDB-free with **~no portal change** (§3 seam) | Low portal-side; long horizon, not portal-controlled | Best *if* it materializes; not actionable today |

**Why (a) is the default and not a cop-out.** The residual is offline-only (§3), the artifact is
upstream-owned (§1), and the build is enormous (§2). Nothing in the portal forces its removal: the
one thing the cube spike wanted — a TileDB-free serving tier — is reachable without touching the
corpus at all.

**Why (b) is heavier than it looks.** The published H5ADs are *pre-integration* raw inputs. WMG's
`value_filter` + `axis_query` run over the *harmonized* corpus (unioned gene universe, normalized X,
consistent ontology terms). Consuming raw H5ADs means owning that harmonization — i.e. maintaining a
fork of CZI's integration in perpetuity, not a one-time format swap. (Cf. the gene-universe work
noted in the findings doc.)

**Why (c) is the worst trade.** It's (b)'s ongoing ownership burden *plus* the full 512-GiB build,
just to change the on-disk format of an artifact that already works. It aligns thematically with the
Explorer `.cxg`→Zarr direction ("corpus is a per-cell tensor → Zarr fits", per the summary), and
that alignment is real — but it only makes sense pushed **upstream into CZI's builder**, where one
change serves every Census consumer, not forked here for WMG alone.

---

## 5. Recommendation

**Keep consuming the upstream TileDB-SOMA corpus (option a).** It is the lazy-correct choice: the
residual is offline-only, upstream-owned, and datacenter-scale, and a TileDB-free *serving* tier
(the actual goal of the cube spike) needs nothing from the corpus.

If a full-stack TileDB purge is ever mandated, the order of preference is:
1. **(d) upstream non-TileDB SOMA backend** — least portal code, one change serves all consumers;
   pursue as an upstream ask/contribution, not a fork.
2. **(b) source-H5AD bypass** — only if (d) never lands and the mandate is hard; accept owning
   integration.
3. **(c) fork the builder** — avoid; only viable as an upstream contribution to CZI's builder.

This reinforces the spike's standing verdict: **keep TileDB today, buy optionality.** The corpus is
the part of the stack the portal least controls and least needs to change — de-risking it means
knowing the paths above exist, not walking one now.
