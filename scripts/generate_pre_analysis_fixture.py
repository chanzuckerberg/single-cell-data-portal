"""
One-time script to generate the pre-analysis h5ad test fixture.

Usage:
    python scripts/generate_pre_analysis_fixture.py [source.h5ad]

If no source path is given the script downloads the standard functional-test
example from Dropbox.  The output is written to:

    tests/functional/backend/fixtures/example_pre_analysis.h5ad

Upload that file to Dropbox (dl=1 link) and set PRE_ANALYSIS_DATASET_URI in
tests/functional/backend/constants.py.
"""

import sys
import urllib.request
from pathlib import Path

import anndata as ad

DROPBOX_DIRECT = (
    "https://www.dropbox.com/scl/fi/l18zfx8i45j90xejip6qk/example_valid.h5ad"
    "?rlkey=3p9cnwc0i0ysnfgk8bqy57fev&st=msygclrx&dl=1"
)

OUTPUT = Path(__file__).parent.parent / "tests/functional/backend/fixtures/example_pre_analysis.h5ad"


def main():
    if len(sys.argv) > 1:
        source = Path(sys.argv[1])
        print(f"Reading {source}")
    else:
        source = Path("/tmp/example_valid.h5ad")
        print(f"Downloading example_valid.h5ad → {source}")
        urllib.request.urlretrieve(DROPBOX_DIRECT, source)

    adata = ad.read_h5ad(source)

    # Strip fields that are not present in pre-analysis datasets.
    adata.obsm = {}
    if "cell_type_ontology_term_id" in adata.obs.columns:
        adata.obs = adata.obs.drop(columns=["cell_type_ontology_term_id"])
    adata.uns.pop("default_embedding", None)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUTPUT)
    print(f"Written → {OUTPUT}")


if __name__ == "__main__":
    main()
