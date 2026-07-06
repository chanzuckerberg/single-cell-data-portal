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

> **Note on ownership.** All three repos in this spike — `single-cell-data-portal`,
> `cellxgene-census`, `single-cell-explorer` — are **CZI** (`chanzuckerberg` org). The boundaries in
> this doc are **component / repo / team** boundaries *inside* CZI, not organizational ones. So a
> corpus change is a **cross-team CZI effort**, not an external dependency or "someone else's call."

- The **cube** is WMG's own artifact: 7 precomputed TileDB arrays of *aggregate* stats
  (`nnz`/`sum`/`sqsum`, cell counts) keyed by gene/tissue/organism/cell-type/etc. Built in the
  **data-portal repo** (`backend/wmg/pipeline/`). The cube spike is entirely about this.
- The **corpus** is the **Census**: the harmonized, integrated per-cell
  expression matrix for all of CELLxGENE, published as a TileDB-SOMA object. Built in a
  **different CZI repo (a separate team)** — `../cellxgene-census/tools/cellxgene_census_builder/`.
  The data-portal repo only *reads* it.

The data-portal builds **no corpus of its own.** WMG opens the Census read-only and aggregates
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

So "moving the Census corpus off TileDB" is **not a change in the data-portal repo** — it's a change
to a datacenter-scale weekly build in another CZI repo (the `cellxgene-census` team) *and* to the
public API that serves it (§3). A cross-team CZI effort. That framing bounds every option in §4.

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

## 3. Where TileDB sits — the portal's read *and* the public API

Two things about the corpus's TileDB dependency, one narrow and one wide.

**Narrow — inside the portal, the corpus read is offline-only.**

- **Serving** (the WMG/DE API tier) does **not** touch the corpus. It reads only the cube. So a
  TileDB-free *serving* tier is already achievable via the chDB cube — no corpus change needed.
- **The corpus read is offline, pipeline-only** — it happens once per snapshot build in the WMG
  pipeline job (`open_soma` above), not on any portal request path.

So *for the portal* the residual is a **build-image dependency**, not a serving one: `tiledbsoma`
must be installed in the offline cube-build job because that's what reads the Census.

**Wide — the TileDB-SOMA object is itself a live public API, not a portal-private file.** The same
corpus the WMG pipeline reads offline is *also* served — live, lazily, over anonymous S3 — as the
public **cellxgene-census reader API**
([docs](https://chanzuckerberg.github.io/cellxgene-census/)). This is a whole consumption surface
*beyond* the portal, and it reads TileDB directly on every call:

- **Python `cellxgene_census` + R `cellxgene.census`.** `open_soma()` opens the corpus straight from
  `s3://cellxgene-data-public/cell-census/<tag>/soma/` via TileDB VFS
  (`vfs.s3.no_sign_request=true`, us-west-2), resolving the build through a public
  `release.json`/`mirrors.json`, and returns a **live `tiledbsoma.Collection`** — no local copy.
  On top: `get_anndata` / `get_seurat` / `get_single_cell_experiment`, `get_obs`/`get_var`,
  `get_presence_matrix`, plus an `experimental` layer (embeddings, PyTorch/Geneformer/HuggingFace ML
  loaders).
- **Broad external userbase, format-as-contract.** R+Python near-parity, ~24 tutorial notebooks, a
  community embeddings-contribution ecosystem (a separate `cellxgene-contrib-public` bucket +
  `contributions.json`), and a versioning model (`stable` LTS / weekly `latest` / `V#`) explicitly
  designed so users can **pin a build for reproducibility**. The reader's `tiledbsoma` version is
  pinned to the builder's specifically for read-compatibility — i.e. the on-disk TileDB-SOMA layout
  is treated as a **stable public contract**.

So the corpus's TileDB isn't only a data-portal build-image detail — it's the storage engine under a
live public API that **many independent clients read directly**. The `single-cell-data-portal` is
**one consumer among many** — including **external, non-CZI researchers** — and that reader API is
owned by another CZI team (the `cellxgene-census` repo), not the data-portal team.

**Why `open_soma` == TileDB, structurally — and why that cuts two ways.** SOMA is an **API spec**
(stack-agnostic by design), but `tiledbsoma` is currently the **only production implementation**, and
the reader API is coupled to it *specifically* — hard `tiledbsoma` dependency, TileDB VFS S3 read
path, public TileDB config keys (`get_default_soma_context`). So today "read the Census via SOMA" and
"read TileDB" are the same statement. This is why a format change (§4 option c) would break far more
than WMG — it breaks the public reader API for the whole community — and simultaneously why a
non-TileDB **SOMA backend** (§4 option d) is the only clean exit: the SOMA abstraction is the one seam
that could carry the entire ecosystem, not just the portal, onto a new store.

---

## 4. Options

The §1 framing — a datacenter-scale weekly build in another CZI repo, serving a public API — is the
filter. "What it takes" is measured in *which CZI team does the work* and *how much*, not in latency
(there's nothing to benchmark).

| # | Option | What it takes | Cost | Verdict |
|---|---|---|---|---|
| a | **Leave it — consume the TileDB-SOMA corpus as-is** | Nothing. `open_soma` stays; TileDB lives only in the offline build image | ~0 | **Recommended** |
| b | **Bypass the corpus — build the cube from source H5ADs** | Ingest the pre-integration H5ADs the Census already publishes (`.../<tag>/h5ads/`), then re-implement the cross-dataset **integration + gene-universe harmonization** to reproduce today's `axis_query` semantics | Large, **data-portal-owned, ongoing** (tracks Census schema/integration changes forever) | Fallback only if forced |
| c | **Change the builder to emit Zarr/Parquet** (obs/var → Parquet, X → sparse Zarr) | In `cellxgene-census`, replace the SOMA/TileDB writer layer in `build_soma/`; still run the 512-GiB weekly build; WMG re-aggregates from the new format | **Largest.** Reworking a datacenter-scale build + breaking the on-disk format | Only done in the `cellxgene-census` repo itself, never a data-portal fork |
| d | **Add a non-TileDB SOMA backend** | If an Arrow/Zarr-backed SOMA implementation lands in the Census stack, `open_soma` becomes TileDB-free with **small, bounded data-portal change** — swap deps + re-point ~4 census-read files off direct `tiledbsoma` symbols (§6); does *not* alone remove TileDB from the data-portal (the WMG cube is separately TileDB) | Low data-portal-side; long horizon, needs the `cellxgene-census` team | Best *if* it materializes; not actionable today |

**Why (a) is the default and not a cop-out.** The residual is offline-only (§3), the artifact is
built and served by another CZI team (§1), and the build is enormous (§2). Nothing in the data-portal
forces its removal: the one thing the cube spike wanted — a TileDB-free serving tier — is reachable
without touching the corpus at all.

**Why (b) is heavier than it looks.** The published H5ADs are *pre-integration* raw inputs. WMG's
`value_filter` + `axis_query` run over the *harmonized* corpus (unioned gene universe, normalized X,
consistent ontology terms). Consuming raw H5ADs means owning that harmonization in the data-portal —
i.e. maintaining a parallel copy of the Census integration in perpetuity, not a one-time format swap.
(Cf. the gene-universe work noted in the findings doc.)

**Why (c) is the worst trade.** It's (b)'s ongoing ownership burden *plus* the full 512-GiB build,
just to change the on-disk format of an artifact that already works — and, per §3, that artifact is a
**public API contract**. A data-portal-side fork that emits Zarr/Parquet either forks the public
`cellxgene_census`/`cellxgene.census` reader too (breaking every external Python/R/embeddings/ML
client) or diverges silently from the corpus the rest of the world reads. It aligns thematically with
the Explorer `.cxg`→Zarr direction ("corpus is a per-cell tensor → Zarr fits", per the summary), and
that alignment is real — but it only makes sense done **in the `cellxgene-census` repo itself (builder
+ reader)**, as a SOMA-backend change (option d), where one change serves every Census consumer. Never
a data-portal fork.

---

## 5. Recommendation

**Keep consuming the TileDB-SOMA corpus (option a).** It is the lazy-correct choice: for the
data-portal the residual is offline-only, and a TileDB-free *serving* tier (the actual goal of the
cube spike) needs nothing from the corpus. And more decisively — the corpus format isn't the
data-portal team's to change alone: it's a **datacenter-scale build in another CZI repo** *and* a
**live public API contract** (§3) read directly by a large external Python/R/ML ecosystem. Any format
move ripples across all of that, not just WMG.

If a full-stack TileDB purge is ever mandated, the order of preference is:
1. **(d) non-TileDB SOMA backend in the Census stack** — least data-portal code, and the *only* path
   that carries the whole ecosystem (data-portal + public reader API) at once via the SOMA seam. Drive
   it with the `cellxgene-census` team, not as a fork.
2. **(b) source-H5AD bypass** — only if (d) never lands and the mandate is hard; accept owning
   integration. Data-portal-local, so it sidesteps (but doesn't help) the public API.
3. **(c) change the builder** — avoid; it either breaks or forks the public reader API. Only viable as
   a `cellxgene-census` builder+reader change, i.e. it collapses into (d).

This reinforces the spike's standing verdict: **keep TileDB today, buy optionality.** The corpus is
the part of the stack the data-portal team least controls and least needs to change — de-risking it
means knowing the paths above exist, not walking one now.

---

## 6. If forced — the high-level scope

The first move is to split what "move the corpus off TileDB" actually means, because two very
different mandates hide in that sentence.

### Mandate A — "drop the TileDB C++ dependency from *our* deployment"

Then you don't touch the corpus at all. It stays TileDB-SOMA; only how the *data-portal* gets its
data changes:
- Serve the cube from **chDB** (already de-risked — see [`WMG_CUBE_TILEDB_EXIT_REWORK.md`](WMG_CUBE_TILEDB_EXIT_REWORK.md)),
  removing TileDB from the *serving* image; and
- for the build image, take the **source-H5AD bypass** (option b): ingest the pre-integration H5ADs
  and re-implement the integration / gene-harmonization in the data-portal.

**Scope:** data-portal-owned, bounded, no cross-team coordination — but the data-portal now owns
cross-dataset integration in perpetuity. This is the realistic "forced" path if the driver is the
data-portal's dependency hygiene.

### Mandate B — "the corpus format itself must not be TileDB anywhere"

This is the big one, and it is **fundamentally a `cellxgene-census` effort, not a data-portal one**
(a different CZI team/repo) — the corpus is both a datacenter-scale build *and* the public reader API
(§3). The load-bearing decision is the **target backend**, which forces a sub-fork:

- **B1 — do it through SOMA (strongly preferred).** SOMA is backend-agnostic; if a non-TileDB SOMA
  backend (Zarr/Arrow-native) exists, `open_soma()` keeps working and the whole ecosystem moves at
  once. The catch: **no production non-TileDB SOMA backend ships today**, so this likely means
  *building or adopting one* — writing a storage engine (SparseNDArray, DataFrame, Collection,
  `AxisQuery` predicate pushdown, lazy remote S3 reads), not doing a migration. **This item dominates
  cost and risk.**
- **B2 — abandon SOMA, bespoke layout** (obs/var → Parquet, X → sparse Zarr). Cheaper to write but
  loses the reader API entirely and must rebuild it (Python + R). Fragments the ecosystem — mentioned
  only to be dismissed.

Assuming **B1**, the work is four tracks plus two closers:

1. **Backend (dominant cost).** A non-TileDB SOMA implementation at production quality: SOMA API over
   Zarr/Arrow, anonymous S3 remote reads, and critically **predicate pushdown on `obs`** (census reads
   filter on `is_primary_data` / `nnz` / `cell_type` — exactly TileDB-SOMA's strength) plus efficient
   sparse `X` slicing. Where the effort lives.
2. **Write side — the builder.** Extend `cellxgene_census_builder`'s writer layer (`build_soma/`,
   `globals.py` platform config) to emit the new backend. The 5-step pipeline and the 512-GiB build
   stay; you swap only the SOMA-write target. With a real SOMA backend this is a re-point, not a
   rewrite.
3. **Read side — the public API.** The reader is coupled to `tiledbsoma` *specifically*
   (`get_default_soma_context` → `SOMATileDBContext`, `vfs.s3.*` keys, pinned `tiledbsoma`).
   Generalize those to be backend-agnostic; `get_anndata`/`get_seurat`/`get_presence_matrix` +
   embeddings/ML loaders ride on SOMA unchanged *if* B1 holds.
4. **Contract & ecosystem migration.** The corpus is a versioned public reproducibility contract
   (`stable`/`latest`/`V#`) — no hard cut. Expect a **dual-publish deprecation window** (both formats
   for N releases), migrating the embeddings-contribution bucket + `contributions.json` +
   `census_contrib` tooling, R/Python parity, and comms for users pinning builds.
5. **Performance re-validation.** The same "does the layout meet read latency" gate the cube spike ran,
   but for a *filtered per-cell tensor* workload (obs-filter + sparse X slice over S3), not an OLAP
   aggregate. Zarr fit the Explorer `.cxg` tensor, but that was dense per-dataset reads; census-scale
   filtered SOMA queries are a different test — the real risk to prove.
6. **Portal payoff — small, but not nil, and not sufficient alone.** Because the portal codes against
   the SOMA abstraction, the code churn is bounded to the **census-read seam**, ~4 files
   (`expression_summary_and_cell_counts.py`, `dataset_metadata.py`, `cell_counts.py`,
   `expression_summary.py`) — swap the `cellxgene-census`/`tiledbsoma` deps for the new SOMA-backend
   package and re-point the direct `tiledbsoma` / `ExperimentAxisQuery` / `soma.AxisQuery` references
   to the generalized SOMA namespace (`somacore`), then re-validate cube parity + read performance.
   **Crucially, this does NOT by itself remove TileDB from the data-portal.** The data-portal has a
   *second*, independent TileDB usage — its own WMG cube (raw `import tiledb` across ~20 files:
   `cube_schema*.py`, `snapshot.py`, `data/tiledb.py`, most `pipeline/*.py`), which a corpus format
   change leaves untouched. Dropping the TileDB **C++ dependency entirely** requires *both* this
   corpus change *and* the separate chDB cube exit ([`WMG_CUBE_TILEDB_EXIT_REWORK.md`](WMG_CUBE_TILEDB_EXIT_REWORK.md));
   either alone leaves the C++ core in the image via the other. B1 is still the only clean exit *for
   the corpus read* — but "clean" means small/bounded data-portal churn, not zero.

### Bottom line

| | Mandate A (data-portal dep hygiene) | Mandate B (corpus format) |
|---|---|---|
| Lead CZI team | data-portal | **`cellxgene-census`** (with data-portal); cross-team, not a data-portal-only call |
| Approach | chDB cube + source-H5AD bypass | non-TileDB SOMA backend (B1) |
| Dominant cost | Owning integration forever | **Building a production SOMA backend** |
| Data-portal code if done | Moderate | Small, bounded (not zero — §6.6) |

If forced tomorrow and the driver is the data-portal's dependency, do **Mandate A**. If it's a true
format mandate, the honest scope is "stand up a non-TileDB SOMA backend in the Census stack" — a
storage-engine project measured in quarters and shared across the community, not a data-portal
migration. It's one CZI program spanning the mapped components; the data-portal team drives it with
the `cellxgene-census` team rather than forking.
