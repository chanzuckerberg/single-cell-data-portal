"""SPIKE: concurrency + memory for chDB vs TileDB on the real WMG cube.

    PYTHONPATH=. venv-311/bin/python scripts/wmg_chdb_spike_concurrency.py <snapshot_dir> <chdb_dir>

The realcube spike measured ONE query at a time in ONE process. Prod is 5 gunicorn/gevent workers in a
16 GB Fargate container. This closes the two unknowns that single-query benchmarking can't see:

  1. MEMORY: each worker holds its own embedded engine + caches, so RAM multiplies by worker count.
     Does 5x the chDB engine fit 16 GB where TileDB's configured tile cache did? (peak aggregate RSS)
  2. LATENCY UNDER LOAD: a native chDB query can't yield to gevent, so concurrent requests to a worker
     queue. Does p95/p99 hold up, or does the great single-query median degrade under offered load?

Models prod with a PROCESS POOL: N workers each open the backend once and pull query jobs off a shared
queue (process-level parallelism = worker count, matching gunicorn). Reports throughput + p50/p95/p99
and peak aggregate RSS, for both backends, sweeping N. Also confirms native calls serialize under
gevent (the mechanism behind risk #2).

Scope: the interactive WMG query path (expression_summary default/indexed/secondary + cell_counts) --
the <10s-budget endpoint where concurrency matters. Diffexp (DE, heavier/less-interactive) is out.
"""

import multiprocessing as mp
import os
import random
import shutil
import sys
import time

import psutil

from backend.common.census_cube.data.criteria import CensusCubeQueryCriteria
from backend.common.census_cube.data.query import CensusCubeQuery, CensusCubeQueryParams
from backend.common.census_cube.data.query_chdb import ChdbCensusCubeQuery
from backend.wmg.api.config import (
    READER_CENSUS_CUBE_CUBE_QUERY_VALID_ATTRIBUTES,
    READER_CENSUS_CUBE_CUBE_QUERY_VALID_DIMENSIONS,
)
from scripts.wmg_chdb_spike_realcube import open_chdb, open_tiledb

TILEDB_WORKERS = [1, 5, 10]  # 5 = prod; 1 = baseline; 10 = overload
CHDB_WORKERS = [1, 3]  # chDB locks its dir -> one DB COPY per worker; capped by local disk (17GB each)
N_JOBS = 2000
MEM_BUDGET_MB = 16000  # prod Fargate task


def _params():
    return CensusCubeQueryParams(
        cube_query_valid_attrs=READER_CENSUS_CUBE_CUBE_QUERY_VALID_ATTRIBUTES,
        cube_query_valid_dims=READER_CENSUS_CUBE_CUBE_QUERY_VALID_DIMENSIONS,
    )


def build_jobs(n, genes, tissues, cell_types, seed=0):
    """A repeatable mix of the four interactive query shapes, randomized filter values per job."""
    rng = random.Random(seed)
    kinds = ["default", "indexed", "secondary", "cell_counts"]
    jobs = []
    for i in range(n):
        kind = kinds[i % len(kinds)]
        job = {"kind": kind, "genes": rng.sample(genes, rng.randint(1, min(4, len(genes))))}
        if kind in ("indexed", "cell_counts"):
            job["tissues"] = rng.sample(tissues, rng.randint(1, min(3, len(tissues))))
        if kind == "secondary":
            job["cell_types"] = rng.sample(cell_types, rng.randint(1, min(3, len(cell_types))))
        jobs.append(job)
    return jobs


def run_job(q, organism, job):
    kind = job["kind"]
    if kind == "default":
        return q.expression_summary_default(
            CensusCubeQueryCriteria(organism_ontology_term_id=organism, gene_ontology_term_ids=job["genes"])
        )
    if kind == "indexed":
        return q.expression_summary(
            CensusCubeQueryCriteria(
                organism_ontology_term_id=organism,
                gene_ontology_term_ids=job["genes"],
                tissue_ontology_term_ids=job["tissues"],
            )
        )
    if kind == "secondary":
        return q.expression_summary(
            CensusCubeQueryCriteria(
                organism_ontology_term_id=organism,
                gene_ontology_term_ids=job["genes"],
                cell_type_ontology_term_ids=job["cell_types"],
            )
        )
    return q.cell_counts(
        CensusCubeQueryCriteria(
            organism_ontology_term_id=organism,
            gene_ontology_term_ids=job["genes"],
            tissue_ontology_term_ids=job["tissues"],
        )
    )


def worker(backend, open_arg, organism, job_q, res_q):
    """open_arg is this worker's own chDB dir (chDB locks per-dir) or the shared TileDB snapshot dir."""
    proc = psutil.Process()
    try:
        if backend == "chdb":
            _, ns = open_chdb(open_arg)
            q = ChdbCensusCubeQuery(ns, _params())
        else:
            ns = open_tiledb(open_arg)
            q = CensusCubeQuery(ns, _params())
    except Exception as e:  # open failure must still post a result or the parent blocks on get()
        res_q.put({"error": f"open: {repr(e)[:180]}"})
        return

    latencies = []
    while True:
        job = job_q.get()
        if job is None:
            break
        t = time.perf_counter()
        try:
            run_job(q, organism, job)
        except Exception as e:  # a worker that can't open/query the shared dir is itself a finding
            res_q.put({"error": repr(e)[:200]})
            return
        latencies.append((time.perf_counter() - t) * 1000)
    res_q.put({"latencies": latencies, "peak_rss_mb": proc.memory_info().rss / 1e6})


def pct(sorted_vals, p):
    if not sorted_vals:
        return float("nan")
    k = min(len(sorted_vals) - 1, int(round(p / 100 * (len(sorted_vals) - 1))))
    return sorted_vals[k]


def load_test(backend, worker_open_args, organism, jobs):
    n_workers = len(worker_open_args)
    job_q, res_q = mp.Queue(), mp.Queue()
    for j in jobs:
        job_q.put(j)
    for _ in range(n_workers):
        job_q.put(None)
    procs = [
        mp.Process(target=worker, args=(backend, open_arg, organism, job_q, res_q))
        for open_arg in worker_open_args
    ]
    t0 = time.perf_counter()
    for p in procs:
        p.start()
    results = [res_q.get() for _ in procs]
    for p in procs:
        p.join()
    wall = time.perf_counter() - t0

    errs = [r["error"] for r in results if "error" in r]
    if errs:
        return {"error": errs[0]}
    lats = sorted(l for r in results for l in r["latencies"])
    rss = sum(r["peak_rss_mb"] for r in results)
    return {
        "throughput": len(lats) / wall,
        "p50": pct(lats, 50),
        "p95": pct(lats, 95),
        "p99": pct(lats, 99),
        "rss_mb": rss,
    }


def gevent_serialization_check(chdb_dir, organism, genes):
    """Confirm risk #2's mechanism: a native chDB query blocks the gevent hub, so greenlets serialize."""
    from gevent import monkey

    monkey.patch_all(thread=False)
    import gevent

    _, ns = open_chdb(chdb_dir)
    q = ChdbCensusCubeQuery(ns, _params())
    crit = lambda: CensusCubeQueryCriteria(organism_ontology_term_id=organism, gene_ontology_term_ids=genes[:3])  # noqa

    for _ in range(5):  # warm caches so the baseline isn't a cold outlier (the confound in the first run)
        q.expression_summary_default(crit())
    warm = []
    for _ in range(5):
        t = time.perf_counter()
        q.expression_summary_default(crit())
        warm.append((time.perf_counter() - t) * 1000)
    single_ms = sorted(warm)[len(warm) // 2]

    G = 10
    t = time.perf_counter()
    gevent.joinall([gevent.spawn(lambda: q.expression_summary_default(crit())) for _ in range(G)])
    wall_ms = (time.perf_counter() - t) * 1000
    factor = wall_ms / single_ms if single_ms else float("nan")
    serialized = factor > G * 0.6
    print(f"\ngevent: warm 1 query {single_ms:.1f}ms; {G} concurrent greenlets {wall_ms:.1f}ms (={factor:.1f}x)"
          f" -> {'SERIALIZED: native chDB call blocks the hub (in-process serving does not overlap)' if serialized else 'overlapped'}")


def main(snapshot_dir, chdb_dir):
    # filter-value pool: pull real ids once so both backends run the identical workload
    sess, ns = open_chdb(chdb_dir)
    organism = sess.query("SELECT organism_ontology_term_id FROM wmg.cell_counts LIMIT 1", "DataFrame").iloc[0, 0]
    org_f = f"organism_ontology_term_id = '{organism}'"

    def distinct(table, col, n):
        return (
            sess.query(f"SELECT DISTINCT {col} FROM wmg.{table} WHERE {org_f} LIMIT {n}", "DataFrame")[col]
            .dropna()
            .tolist()
        )

    genes = distinct("expression_summary", "gene_ontology_term_id", 40)
    tissues = distinct("cell_counts", "tissue_ontology_term_id", 15)
    cell_types = distinct("cell_counts", "cell_type_ontology_term_id", 15)
    sess.close()
    jobs = build_jobs(N_JOBS, genes, tissues, cell_types)

    # chDB locks its data dir, so N workers need N separate DB copies (the "5 workers share one synced
    # snapshot" model does NOT work). Pre-make (max-1) copies; disk caps how many workers we can test.
    max_chdb = max(CHDB_WORKERS)
    chdb_dirs = [chdb_dir]
    for i in range(1, max_chdb):
        dst = f"{chdb_dir}_copy{i}"
        if not os.path.isdir(dst):
            print(f"copying chDB DB -> {dst} (chDB can't share a dir across processes) ...")
            shutil.copytree(chdb_dir, dst)
        chdb_dirs.append(dst)

    print(f"\norganism={organism}  {N_JOBS} jobs  mem budget={MEM_BUDGET_MB}MB")
    print(f"{'backend':<8} {'workers':>7} {'thru q/s':>10} {'p50 ms':>8} {'p95 ms':>8} {'p99 ms':>8} {'peak RSS MB':>12}")
    print("-" * 70)
    plan = [("tiledb", n, [snapshot_dir] * n) for n in TILEDB_WORKERS]
    plan += [("chdb", n, chdb_dirs[:n]) for n in CHDB_WORKERS]
    for backend, n, open_args in plan:
        r = load_test(backend, open_args, organism, jobs)
        if "error" in r:
            print(f"{backend:<8} {n:>7}  ERROR: {r['error']}")
            continue
        budget = "" if r["rss_mb"] <= MEM_BUDGET_MB else "  OVER 16GB!"
        print(
            f"{backend:<8} {n:>7} {r['throughput']:>10.0f} {r['p50']:>8.1f} {r['p95']:>8.1f} "
            f"{r['p99']:>8.1f} {r['rss_mb']:>12.0f}{budget}"
        )

    gevent_serialization_check(chdb_dir, organism, genes)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    mp.set_start_method("spawn")
    main(sys.argv[1], sys.argv[2])
