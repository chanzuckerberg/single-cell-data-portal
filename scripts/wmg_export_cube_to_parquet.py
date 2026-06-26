"""SPIKE helper: export TileDB WMG cubes to Parquet (streaming, memory-bounded).

    PYTHONPATH=. python scripts/wmg_export_cube_to_parquet.py <snapshot_dir> <parquet_out_dir>

Streams each cube via TileDB's incomplete-query iterator so the ~4GB expression_summary
never fully materializes in memory. Dims (the TileDB index) are reset to columns so DuckDB
sees a flat table.
"""

import os
import sys
import time

import pyarrow as pa
import pyarrow.parquet as pq
import tiledb

from backend.common.census_cube.data.tiledb import create_ctx

CUBES = ["expression_summary", "expression_summary_default", "cell_counts", "marker_genes"]


def export_cube(cube_dir: str, out_path: str) -> int:
    ctx = create_ctx()
    rows = 0
    writer = None
    with tiledb.open(cube_dir, ctx=ctx) as cube:
        # incomplete iterator -> bounded batches (size set by py.init_buffer_bytes in the ctx)
        for batch in cube.query(return_incomplete=True, use_arrow=False).df[:]:
            if batch.empty:
                continue
            df = batch.reset_index()
            # TileDB ascii dims/attrs come back as bytes; decode so Parquet stores UTF-8 strings
            # (else DuckDB sees BLOBs and string filters won't match).
            for col in df.columns:
                if df[col].dtype == object:
                    df[col] = df[col].map(lambda v: v.decode() if isinstance(v, (bytes, bytearray)) else v)
            table = pa.Table.from_pandas(df, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(out_path, table.schema, compression="zstd")
            writer.write_table(table)
            rows += len(batch)
    if writer is not None:
        writer.close()
    return rows


def main(snapshot_dir: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)
    for name in CUBES:
        cube_dir = os.path.join(snapshot_dir, name)
        if not os.path.isdir(cube_dir):
            print(f"skip {name} (not found)")
            continue
        out_path = os.path.join(out_dir, f"{name}.parquet")
        t = time.perf_counter()
        n = export_cube(cube_dir, out_path)
        mb = os.path.getsize(out_path) / 1024 / 1024
        print(f"{name}: {n:,} rows -> {out_path} ({mb:.0f} MB) in {time.perf_counter()-t:.1f}s")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
