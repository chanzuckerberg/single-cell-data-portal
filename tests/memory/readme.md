Befor running, make sure the containers have INSTALL_DEV=true. This will install the necessary packages for memory profiling. Also install the tools locally to make analysis easier.

TODO: run all anylsis in the container.

## tested datasets:

These are the dataset that were used. They are copied to the test fixture folder here `tests/unit/backend/fixtures`.

- 10x https://datasets.cellxgene.cziscience.com/f50deffa-43ae-4f12-85ed-33e45040a1fa.h5ad
- slide seq https://datasets.cellxgene.cziscience.com/9c889e34-c6b7-4493-affa-be30ca85a299.h5ad
- visium https://datasets.cellxgene.cziscience.com/f9dc9230-cd78-4cc4-8b05-b47341e6c977.h5ad

## for each dataset run the following tests:

### Run memory-profiler and produce a plot:

```bash
docker compose --project-directory ../../ run --rm -w /single-cell-data-portal processing mprof run ./tests/memory/processing/test_process_cxg.py
```

```bash
mprof plot
```

### Run memory-profiler and get a line by line analysis:

```bash
docker compose --project-directory ../../ run --rm -w /single-cell-data-portal processing python -m memory_profiler ./tests/memory/processing/test_process_cxg.py
```

### Run memray and produce a flamegraph:

Before running you will need to comment out @profile from some of the code. This is used by memory-profiler and is not supported by memray.

```bash
docker compose --project-directory ../../ run --rm -w /single-cell-data-portal processing memray run -o memray-test_make_cxg.bin -m tests.memory.processing.test_process_cxg
```

```bash
memray flamegraph memray-test_make_cxg.bin
```
