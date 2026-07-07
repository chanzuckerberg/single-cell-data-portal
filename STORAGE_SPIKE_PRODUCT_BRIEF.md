# Moving off TileDB — options & path forward (findings for product)

*An analysis of whether and how the single-cell platform can move off TileDB, its foundational storage
technology — with a proven path and rough cost per component. Research spike, mid-2026, ~4 min read.*

## TL;DR

- **Moving off TileDB is feasible for all three parts of the platform** — we validated a concrete
  replacement path for each.
- **It's one effort across three CZI components** (three repos). The Census corpus is the heaviest and
  is shared with the outside research community, so its part is a cross-team decision.
- **The "how" is de-risked** — proven paths and rough cost per component — so a decision on *whether*
  and *when* to move can be made on facts rather than guesswork.

*(Whether to prioritize a move, and when, is a product/leadership call. This brief gives the options
and costs to inform it; it doesn't make that call.)*

## Why consider moving off TileDB

We scoped the move independent of the driver — the paths below hold regardless of *why*. That said, the
reasons an org would take this on:

- **It's a heavy, specialized dependency.** TileDB is a native (C++) storage engine; removing it
  simplifies our build and deployment.
- **The ecosystem is moving to open, cloud-native storage.** Formats like Zarr / Arrow / SQL are
  becoming the norm across single-cell tooling and ML pipelines; aligning improves interoperability.
- **Strategic optionality.** Not being locked to one specialized vendor format for our core data.

## The system at a glance

![System map: three CZI components store data in TileDB today; users reach them via the WMG/DE API and the Explorer app](storage_spike_system_map.png)

*Three components store their data in TileDB today (▣). All three live in CZI repos. Users reach the
platform through the WMG / Differential-Expression API and the Explorer app — and the WMG cube is
built from the Census corpus, which is also a public dataset used by researchers worldwide.*

## What we analyzed — three parts, three paths

| Component | Where it lives (CZI) | Replacement path | Status | Effort |
|---|---|---|---|---|
| **Explorer matrix** — data behind the cell viewer | single-cell-explorer | Open format (Zarr) | ✅ Proven | Low–Medium |
| **WMG cube** — powers Where's My Gene + Diff-Expr | single-cell-data-portal | Validated engine swap | ✅ Validated, not built | Medium |
| **Census corpus** — the source dataset the cube is built from | cellxgene-census | Open storage backend | 🟡 Feasible; heaviest & shared | Large (cross-team + public) |

- **Explorer matrix** — the per-dataset data behind the cell viewer. A modern, open-format replacement
  is proven.
- **WMG cube** — powers "Where's My Gene" and Differential Expression. We validated a replacement
  approach that matches its current speed; it is scoped but not yet built.
- **Census corpus** — the large, harmonized dataset the cube is built from. Produced by a different CZI
  team and *also a public product* used by researchers worldwide. We're one consumer of it, so changing
  its format is a shared, cross-team effort with the widest impact.

## The path forward

Treat a move as one program with three workstreams, in order of how independently we can do them:

- **Explorer + WMG cube** — ours to execute; medium effort each; the paths are proven, so these can
  proceed on our own timeline.
- **Census corpus** — the heavy, shared workstream, with two flavors:
  - *Remove the dependency from our own services only:* doable on our side without changing the Census,
    at the cost of taking on some data-preparation work ourselves.
  - *Move the Census itself to an open backend:* a large, shared effort with the cellxgene-census team
    (it's a huge weekly build and a public API), best done as an open storage backend that keeps every
    existing user working.
- **Removing TileDB *entirely* requires all three** — no single change does it.

## Decision & next steps

- **The go/no-go and timing are a product/leadership decision.** This analysis provides the options,
  feasibility, and rough cost to make it.
- **If we proceed:** start with the components we own and have proven (Explorer, WMG cube), and open a
  conversation with the cellxgene-census team on an open corpus backend (the long-lead item).
- **If we hold:** the paths are documented and ready, so the work becomes execution rather than
  research whenever a driver makes it a priority.

## Go deeper (technical detail)

Full scoping and evidence:
[Technical summary](STORAGE_SPIKE_TECHNICAL_SUMMARY.md) ·
[WMG cube exit plan](WMG_CUBE_TILEDB_EXIT_REWORK.md) ·
[Engine architecture](CLICKHOUSE_VS_TILEDB_ARCHITECTURE.md) ·
[Census corpus scope](CENSUS_CORPUS_GENERATION_SCOPE.md)

*For the curious: the validated cube replacement is embedded ClickHouse (served via a ClickHouse
sidecar); the open corpus standard is SOMA over Zarr/Arrow. Specifics are in the linked docs.*
