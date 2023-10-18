Scripts herein are used to verify compatibility between the Python and Typescript implementations of the diffexpdu codec.

This directory includes:

- codecTest.py -- wraps the Python implementation
- codecTest.ts -- wraps the Typescript implementation
- genTestSpec.py -- generates test cases from H5AD files.

To use, you must pip install a few more packages than currently used by explorer backend:

- jsonlines

and for Typescript, have npm modules installed:

- @types/node
- commander
- ts-node (or some other means of running TS)

Examples (all assuming you are running from this directory):

1. Generate a test spec from one H5AD:

   ```
   $ PYTHONPATH=../.. python3 genTestSpec.py --out spec.jsonl.gz datafile.h5ad
   Spec written to spec.jsonl.gz.
   ```

2. Encode with Python implementation, and verify that the same implementation can read it back.

   ```
   $ PYTHONPATH=../.. python3 codecTest.py encode spec.jsonl.gz -b py_outdir
   All encoded.
   $ PYTHONPATH=../.. python3 codecTest.py verify spec.jsonl.gz -b py_outdir
   All verified.
   ```

3. Encode & verify with Typescript implementation:

```
$ npx ts-node codecTest.ts encode spec.jsonl.gz -b ts_outdir
All encoded.
$ npx ts-node codecTest.ts verify spec.jsonl.gz -b ts_outdir
All verified.
```

4. Cross-implementation verification, to ensure Python can read Typescript and vice-versa
   (note change in directory names).

```
$ npx ts-node codecTest.ts verify spec.jsonl.gz -b py_outdir
All verified.
$ PYTHONPATH=../.. python3 codecTest.py verify spec.jsonl.gz -b ts_outdir
All verified.
```
