import argparse
import gzip
import os
import os.path
import sys

import jsonlines
import numpy as np

from server.common.diffexpdu import DiffExArguments


def main():
    parser = argparse.ArgumentParser(description="Generate test specification from H5AD(s).")
    parser.add_argument("--verbose", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--limit", type=int, default=None, help=argparse.SUPPRESS)

    subparsers = parser.add_subparsers(dest="command")

    sp = subparsers.add_parser("encode", help="Create block encoded files from test specification.")
    sp.add_argument("tests", type=str, help="Test specification file name.")
    sp.add_argument("--bindir", "-b", type=str, help="Encoded file directory (output location)")

    sp = subparsers.add_parser("verify", help="Verify block encoded files from test specification.")
    sp.add_argument("--bindir", "-b", type=str, help="Encoded file directory (input location).")
    sp.add_argument("tests", type=str, help="Test specification file name.")

    args = parser.parse_args()

    if args.command == "encode":
        return encode_tests(args)
    if args.command == "verify":
        return verify_tests(args)

    print("Unknown command")
    parser.print_help()
    return 1


def encode_tests(args):
    bindir = args.bindir
    if not os.path.exists(bindir):
        os.mkdir(bindir)
    elif not os.path.isdir(bindir):
        print("Binary (output) directory name already exists and is not a directory")
        return 1

    with gzip.open(args.tests, mode="rb") if args.tests.endswith(".gz") else open(
        args.tests, mode="rb"
    ) as fp, jsonlines.Reader(fp) as tests:
        limit = 0
        for test in tests:
            test_id = test["test_id"]
            if args.verbose:
                print(f"{test_id}")
            de_args = DiffExArguments(
                mode=DiffExArguments.DiffExMode.TopN,
                params=DiffExArguments.TopNParams(N=50),
                set1=test["set1"]["postings_list"],
                set2=test["set2"]["postings_list"],
            )
            encoded_postings = de_args.pack()
            with open(os.path.join(bindir, f"{test_id}.bin"), mode="wb") as bin:
                bin.write(encoded_postings)

            limit += 1
            if args.limit and limit > args.limit:
                break

    print("All encoded.")
    return 0


def verify_tests(args):
    bindir = args.bindir
    if not os.path.isdir(bindir):
        print("Binary file directory not found.")
        return 1

    with gzip.open(args.tests, mode="rb") if args.tests.endswith(".gz") else open(
        args.tests, mode="rb"
    ) as fp, jsonlines.Reader(fp) as tests:
        limit = 0
        for test in tests:
            test_id = test["test_id"]
            if args.verbose:
                print(f"{test_id}")

            with open(os.path.join(bindir, f"{test_id}.bin"), mode="rb") as bin:
                buf = bin.read()

            de = DiffExArguments.unpack_from(buf)
            assert de.mode == DiffExArguments.DiffExMode.TopN
            assert de.params.N == 50
            assert np.array_equal(de.set1, test["set1"]["postings_list"])
            assert np.array_equal(de.set2, test["set2"]["postings_list"])

            limit += 1
            if args.limit and limit > args.limit:
                break

    print("All verified.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
