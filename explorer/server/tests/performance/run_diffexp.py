import argparse
import gc
import sys
import time

import numpy as np

from server.common.config.app_config import AppConfig
from server.dataset.matrix_loader import DataLoader


def main():
    parser = argparse.ArgumentParser("A command to test diffexp")
    parser.add_argument("dataset", help="name of a dataset to load")
    parser.add_argument("-n", "--num", nargs=2, type=int, help="number of rows to draw for set A and B")
    parser.add_argument("-v", "--var", nargs=2, help="obs variable:value for set A and B")
    parser.add_argument("-f", "--fraction", nargs=2, type=float, help="fraction (0, 1] of obs to drat for set A and B")

    parser.add_argument("-t", "--trials", default=1, type=int, help="number of trials")
    parser.add_argument("-a", "--alg", choices=("cxg",), default="cxg", help="algorithm to use")
    parser.add_argument("-s", "--show", default=False, action="store_true", help="show the results")
    parser.add_argument(
        "--new-selection", default=False, action="store_true", help="change the selection between each trial"
    )
    parser.add_argument("--seed", default=1, type=int, help="set the random seed")

    args = parser.parse_args()
    if sum([int(bool(arg)) for arg in [args.num, args.var, args.fraction]]) != 1:
        print("Must supply one and only one of --num, --var or --fraction")
        return 1
    if args.var and args.new_selection:
        print("--new-selection and --var are incompatible.")
        return 1

    app_config = AppConfig()
    app_config.update_server_config(app__verbose=True)
    app_config.update_server_config(app__flask_secret_key="howdy")
    # CXG Adaptor config - directly influences diffex performance.
    # Need to ensure these are correctly set in the deployment config.
    app_config.update_server_config(
        adaptor__cxg_adaptor__tiledb_ctx={
            "sm.tile_cache_size": "60129542144",
            "py.init_buffer_bytes": str(512 * 1024**2),
            "vfs.s3.region": "us-west-2",
        }
    )
    app_config.complete_config()

    loader = DataLoader(location=args.dataset, app_config=app_config)
    adaptor = loader.open()
    np.set_printoptions(edgeitems=10, linewidth=180)
    rng = np.random.default_rng(seed=args.seed)

    if args.show:
        adaptor.open_X_array().schema.dump()

    print(f"Dataset shape: {adaptor.get_shape()}")
    print(f"Sparse: {adaptor.open_X_array().schema.sparse}")

    # The _first_ call to compute differential expression will compile (numba) various
    # functions. We want to remove that from the benchmark, as it only happens once
    # on backend startup.
    adaptor.compute_diffexp_ttest(np.array([0]), np.array([1]), selector_lists=True)

    filterA, filterB = draw_cell_sets(args, adaptor, rng)

    times = []
    for _i in range(args.trials):
        gc.collect()
        t1 = time.time()
        if args.alg == "cxg":
            results = adaptor.compute_diffexp_ttest(filterA, filterB, selector_lists=True)
        else:
            raise ValueError(f"Unsupported algo {args.alg}")

        t2 = time.time()
        print("TIME=", t2 - t1)
        times.append(t2 - t1)

        if args.new_selection:
            filterA, filterB = draw_cell_sets(args, adaptor, rng)

    if args.show:
        print(results.get("positive", [])[:3])

    print(f"time mean={np.mean(times)}, std={np.std(times)}, n_obs={len(filterA)+len(filterB)}")

    return 0


def draw_cell_sets(args, adaptor, rng):
    n_obs = adaptor.get_shape()[0]

    if args.num:
        assert sum(args.num) <= n_obs, "--num exceeds n_obs"
        draw = rng.choice(n_obs, size=sum(args.num), replace=False)
        return draw[0 : args.num[0]], draw[args.num[0] :]

    if args.var:
        vname, vval = args.var[0].split(":")
        filterA = get_filter_from_obs(adaptor, vname, vval)
        vname, vval = args.var[1].split(":")
        filterB = get_filter_from_obs(adaptor, vname, vval)
        return filterA, filterB

    if args.fraction:
        assert 0 < sum(args.fraction) <= 1.0, "--fraction not in range (0, 1]"
        num = [max(1, int(n_obs * f)) for f in args.fraction]
        draw = rng.choice(n_obs, size=sum(num), replace=False)
        return draw[0 : num[0]], draw[num[0] :]

    # error
    raise AssertionError("Unexpected cli arguments")


def get_filter_from_obs(adaptor, obsname, obsval):
    attrs = adaptor.get_obs_columns()
    if obsname not in attrs:
        print(f"Unknown obs attr {obsname}: expected on of {attrs}")
        sys.exit(1)
    obsvals = adaptor.query_obs_array(obsname)[:]
    obsval = type(obsvals[0])(obsval)

    vfilter = np.where(obsvals == obsval)[0]
    if len(vfilter) == 0:
        u = np.unique(obsvals)
        print(f"Unknown value in variable {obsname}:{obsval}: expected one of {list(u)}")
        sys.exit(1)

    return vfilter


if __name__ == "__main__":
    sys.exit(main())
