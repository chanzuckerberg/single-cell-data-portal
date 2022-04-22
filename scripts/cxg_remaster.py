import argparse
import time

import numpy as np
import psutil
import tiledb

from backend.corpora.dataset_processing.remaster_cxg import compute

"""
Re-master an existing "X" tiledb array in a CXG
Creates a new array with the same data but a different tdb schema
"""

X_extent = None
X_name = None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cxg", type=str, help="cxg path")
    parser.add_argument(
        "--kind", choices=["auto", "sparse", "dense"], type=str, help="target: auto, sparse or dense", default="auto"
    )
    parser.add_argument("--obs", type=int, help="obs extent", default=256)
    parser.add_argument("--var", type=int, help="var extent", default=2048)
    parser.add_argument("--cell-order", choices=["row", "col"], type=str, default="row")
    parser.add_argument("--tile-order", choices=["row", "col"], type=str, default="col")
    parser.add_argument("--capacity", type=int, default=128000)
    parser.add_argument("--compression", type=int, default=22)
    parser.add_argument("--target-array", type=str, default="X_new")
    parser.add_argument("--source-array", type=str, default="X")
    parser.add_argument("--sparse-threshold", type=float, default=25.0, help="sparse encoding threshold")
    args = parser.parse_args()

    compute(    
        args.cxg,
        args.obs,
        args.var,
        args.source_array,
        args.target_array,
        args.kind,
        args.sparse_threshold,
        args.compression,
        args.cell_order,
        args.tile_order,
        args.capacity,
    )

if __name__ == "__main__":
    main()
