"""SPIKE: the SIDECAR serving model — clickhouse-server + socket clients, vs embedded chDB.

    PYTHONPATH=. venv-311/bin/python scripts/wmg_ch_sidecar_concurrency.py [host] [port]

The embedded chDB concurrency spike (wmg_chdb_spike_concurrency.py) found the fatal issues for
in-process serving: exclusive dir lock, per-worker DB copies, and gevent queries serializing ~15x
(a native call blocks the hub). The proposed fix is a clickhouse-server sidecar the web workers hit
over localhost, so a query becomes a SOCKET CALL that gevent can yield on. This harness proves that:

  A. gevent overlap: G greenlets issuing concurrent queries to the server should finish in ~1x a
     single query (the hub schedules the others while each awaits its socket) -- vs embedded's ~15x.
  B. process-pool load: N client processes (= gunicorn workers) hitting ONE server; p50/p95/p99 +
     throughput, and the server is ONE process with SHARED caches (server RSS sampled separately in
     the orchestration, not multiplied by workers).

Assumes a clickhouse-server already running with the wmg tables loaded (see the orchestration that
starts it). Queries mirror the four interactive shapes from the embedded spike.
"""

import multiprocessing as mp
import random
import sys
import time

DB = "wmg"


def _in(vals):
    return "(" + ", ".join("'" + str(v).replace("'", "''") + "'" for v in vals) + ")"


def build_sql(job, organism):
    kind = job["kind"]
    if kind == "cell_counts":
        clauses = [f"organism_ontology_term_id = '{organism}'", f"tissue_ontology_term_id IN {_in(job['tissues'])}"]
        return f"SELECT * FROM {DB}.cell_counts WHERE {' AND '.join(clauses)}"
    table = "expression_summary_default" if kind == "default" else "expression_summary"
    clauses = [f"organism_ontology_term_id = '{organism}'", f"gene_ontology_term_id IN {_in(job['genes'])}"]
    if kind == "indexed":
        clauses.append(f"tissue_ontology_term_id IN {_in(job['tissues'])}")
    if kind == "secondary":
        clauses.append(f"cell_type_ontology_term_id IN {_in(job['cell_types'])}")
    return f"SELECT * FROM {DB}.{table} WHERE {' AND '.join(clauses)}"


def client(host, port):
    import clickhouse_connect

    return clickhouse_connect.get_client(host=host, port=port)


def build_jobs(n, genes, tissues, cell_types, seed=0):
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


def pct(sorted_vals, p):
    if not sorted_vals:
        return float("nan")
    return sorted_vals[min(len(sorted_vals) - 1, int(round(p / 100 * (len(sorted_vals) - 1))))]


# ---- Part A: gevent overlap (the load-bearing proof) --------------------------------------------

def gevent_overlap(host, port, organism, genes):
    from gevent import monkey

    monkey.patch_all(thread=False)
    import gevent

    job = {"kind": "default", "genes": genes[:3]}
    sql = build_sql(job, organism)

    def one():
        c = client(host, port)  # own connection per greenlet -> N concurrent HTTP requests
        c.query(sql)

    for _ in range(3):  # warm
        one()
    t = time.perf_counter()
    one()
    single_ms = (time.perf_counter() - t) * 1000

    G = 10
    t = time.perf_counter()
    gevent.joinall([gevent.spawn(one) for _ in range(G)])
    wall_ms = (time.perf_counter() - t) * 1000
    factor = wall_ms / single_ms if single_ms else float("nan")
    verdict = "OVERLAPPED: socket call yields to the hub (embedded chDB serialized ~15x here)" if factor < G * 0.6 else "serialized"
    print(f"\ngevent: warm 1 query {single_ms:.1f}ms; {G} concurrent greenlets {wall_ms:.1f}ms (={factor:.1f}x) -> {verdict}")


# ---- Part B: process-pool load (gunicorn workers as clients, one shared server) -----------------

def worker(host, port, organism, job_q, res_q):
    c = client(host, port)
    lat = []
    while True:
        job = job_q.get()
        if job is None:
            break
        t = time.perf_counter()
        c.query(build_sql(job, organism))
        lat.append((time.perf_counter() - t) * 1000)
    res_q.put(lat)


def load_test(host, port, organism, jobs, n_clients):
    job_q, res_q = mp.Queue(), mp.Queue()
    for j in jobs:
        job_q.put(j)
    for _ in range(n_clients):
        job_q.put(None)
    procs = [mp.Process(target=worker, args=(host, port, organism, job_q, res_q)) for _ in range(n_clients)]
    t0 = time.perf_counter()
    for p in procs:
        p.start()
    results = [res_q.get() for _ in procs]
    for p in procs:
        p.join()
    wall = time.perf_counter() - t0
    lats = sorted(l for r in results for l in r)
    return {"throughput": len(lats) / wall, "p50": pct(lats, 50), "p95": pct(lats, 95), "p99": pct(lats, 99)}


def main(host, port):
    c = client(host, port)
    organism = c.query(f"SELECT organism_ontology_term_id FROM {DB}.cell_counts LIMIT 1").result_rows[0][0]
    org_f = f"organism_ontology_term_id = '{organism}'"

    def distinct(table, col, n):
        rows = c.query(f"SELECT DISTINCT {col} FROM {DB}.{table} WHERE {org_f} LIMIT {n}").result_rows
        return [r[0] for r in rows if r[0]]

    genes = distinct("expression_summary", "gene_ontology_term_id", 40)
    tissues = distinct("cell_counts", "tissue_ontology_term_id", 15)
    cell_types = distinct("cell_counts", "cell_type_ontology_term_id", 15)
    jobs = build_jobs(2000, genes, tissues, cell_types)
    print(f"organism={organism}  server={host}:{port}  2000 jobs")

    print(f"\n{'clients':>7} {'thru q/s':>10} {'p50 ms':>8} {'p95 ms':>8} {'p99 ms':>8}")
    print("-" * 46)
    for n in (1, 5, 10):
        r = load_test(host, port, organism, jobs, n)
        print(f"{n:>7} {r['throughput']:>10.0f} {r['p50']:>8.1f} {r['p95']:>8.1f} {r['p99']:>8.1f}")

    gevent_overlap(host, port, organism, genes)


if __name__ == "__main__":
    host = sys.argv[1] if len(sys.argv) > 1 else "localhost"
    port = int(sys.argv[2]) if len(sys.argv) > 2 else 8123
    mp.set_start_method("spawn")
    main(host, port)
