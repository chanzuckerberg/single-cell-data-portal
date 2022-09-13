import argparse

import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.dataset_processing.remaster_cxg_v2 import compute

"""
Re-master an existing "X" tiledb array in a CXG
Creates two new arrays (X and Xc) with the same data but row- and column-optimized schemas
"""

X_extent = None
X_name = None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cxg", type=str, help="cxg path")
    parser.add_argument(
        "--kind", choices=["auto", "sparse", "dense"], type=str, help="target: auto, sparse or dense", default="auto"
    )
    parser.add_argument("--obs-extent", type=int, help="obs extent", default=256)
    parser.add_argument("--var-extent", type=int, help="var extent", default=256)
    parser.add_argument("--cell-order", choices=["row", "col"], type=str, default="row")
    parser.add_argument("--tile-order", choices=["row", "col"], type=str, default="row")
    parser.add_argument("--capacity", type=int, default=1024000)
    parser.add_argument("--compression", type=int, default=22)
    parser.add_argument("--target-array", type=str, default="X_new")
    parser.add_argument("--source-array", type=str, default="X")
    args = parser.parse_args()

    compute(**vars(args))


if __name__ == "__main__":
    main()
