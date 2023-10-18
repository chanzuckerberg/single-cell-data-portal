import argparse
import gc
import gzip
import itertools
import os
import os.path
import sys

import anndata
import jsonlines
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Generate test specification from H5AD(s).")
    parser.add_argument("--verbose", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("h5ad", nargs="+", help="H5AD file.")
    parser.add_argument("--seed", type=int, default=0, help="random number seed")
    parser.add_argument(
        "--out", "-o", type=str, help="test specification output file name. Will append '.gz' if not specified."
    )
    parser.add_argument("--trials", type=int, default=1, help="number of random trials to generate per random test.")
    args = parser.parse_args()

    return generate_tests(args)


def generate_tests(args):
    """Write tests in JSON format to output file."""

    fname = args.out
    if not fname.endswith(".jsonl.gz"):
        if fname.endswith(".jsonl"):
            fname += ".gz"
        else:
            fname += ".jsonl.gz"

    with gzip.open(fname, mode="wb") as fp, jsonlines.Writer(fp) as output:
        rng = np.random.default_rng(args.seed)
        test_id = 0

        for h5ad in args.h5ad:
            if args.verbose:
                print("loading", h5ad)
            ad = anndata.read_h5ad(h5ad, backed=True)

            # Reject the H5AD if it doesn't contain at least some useful obs annotation.
            if not any(cat in ad.obs for cat in ["cell_type", "tissue", "sex"]):
                print(f"Sorry, no useful categories in obs for file {h5ad} -- skipping this file.")
                continue

            n_obs = ad.n_obs
            dataset_id = os.path.splitext(os.path.basename(h5ad))[0]

            for (fracA, fracB) in itertools.product([0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.49], repeat=2):
                n_elem = (int(n_obs * fracA), int(n_obs * fracB))
                if args.verbose:
                    print(f"{dataset_id}: uniform {fracA}, {fracB}")
                for trial in range(args.trials):
                    draw = rng.choice(n_obs, sum(n_elem), replace=False).astype(np.uint32)  # uniform distribution
                    set1 = {"label": str(fracA), "postings_list": draw[0 : n_elem[0]].tolist()}
                    set2 = {"label": str(fracB), "postings_list": draw[n_elem[0] :].tolist()}
                    set1["postings_list"].sort()
                    set2["postings_list"].sort()
                    if args.verbose:
                        print(f"{test_id}: {dataset_id}: uniform, {str(fracA)}, {str(fracB)}")
                    output.write(
                        {
                            "test_id": test_id,
                            "dataset_id": dataset_id,
                            "cat": "uniform",
                            "trial": trial,
                            "set1": set1,
                            "set2": set2,
                        }
                    )
                    test_id += 1

            # real data labels encoded as diffex params.
            for col in ["tissue", "sex", "cell_type"]:
                labels = ad.obs[col].unique().tolist()
                for (labelA, labelB) in itertools.combinations(labels, 2):
                    setA = {
                        "label": labelA,
                        "postings_list": (np.asarray(ad.obs[col] == labelA).nonzero()[0]).astype(np.uint32).tolist(),
                    }
                    setB = {
                        "label": labelB,
                        "postings_list": (np.asarray(ad.obs[col] == labelB).nonzero()[0]).astype(np.uint32).tolist(),
                    }
                    setA["postings_list"].sort()
                    setB["postings_list"].sort()
                    if args.verbose:
                        print(f"{test_id}: {dataset_id}: {col}, {labelA}, {labelB}")
                    output.write(
                        {
                            "test_id": test_id,
                            "dataset_id": dataset_id,
                            "cat": col,
                            "trial": 0,
                            "set1": set1,
                            "set2": set2,
                        }
                    )
                    test_id += 1

            # now use "real" data labels for selection of draws. This draws each label in each category,
            # AND each combination (of 2) in each category.
            for col in ["cell_type", "tissue", "sex"]:
                labels = ad.obs[col].unique().tolist()
                test_labels = list(itertools.combinations(labels, 1)) + list(itertools.combinations(labels, 2))
                for test_label in test_labels:
                    postings_list = np.isin(ad.obs[col], test_label).nonzero()[0]
                    label_description = " & ".join(test_label)
                    postings_list = postings_list.astype(np.uint32)
                    if args.verbose:
                        print(f"{test_id}: {dataset_id}: {col}, {label_description}")
                    output.write(
                        {
                            "test_id": test_id,
                            "dataset_id": dataset_id,
                            "cat": col,
                            "trial": 0,
                            "set1": {"label": label_description, "postings_list": postings_list.tolist()},
                            "set2": {"label": "cell id 0", "postings_list": [0]},
                        }
                    )
                    test_id += 1

            del ad
            gc.collect()

    print(f"Spec written to {fname}.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
