import os

API_URL = {
    "prod": "https://api.cellxgene.cziscience.com",
    "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
    "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
    "test": "https://localhost:5000",
    "rdev": f"https://{os.getenv('STACK_NAME', '')}-backend.rdev.single-cell.czi.technology",
}
AUDIENCE = {
    "prod": "api.cellxgene.cziscience.com",
    "staging": "api.cellxgene.staging.single-cell.czi.technology",
    "test": "api.cellxgene.dev.single-cell.czi.technology",
    "dev": "api.cellxgene.dev.single-cell.czi.technology",
    "rdev": "api.cellxgene.dev.single-cell.czi.technology",
}

DATASET_URI = "https://www.dropbox.com/scl/fi/l18zfx8i45j90xejip6qk/example_valid.h5ad?rlkey=3p9cnwc0i0ysnfgk8bqy57fev&st=msygclrx&dl=0"

VISIUM_DATASET_URI = "https://www.dropbox.com/scl/fi/a871dhnoqzzhkqz0pr3fj/small_visium.h5ad?rlkey=8zv1m79p4lnaqr60nl6hi6yai&st=ibim78id&dl=0"

DATASET_MANIFEST = {"anndata": DATASET_URI}
VISIUM_DATASET_MANIFEST = {"anndata": VISIUM_DATASET_URI}
ATAC_SEQ_MANIFEST = {
    "anndata": "https://www.dropbox.com/scl/fi/1h9d1usgobiqgnin4s6ld/atac.h5ad?rlkey=o3aokwmez6a8tdip5p4vybv8a&st=nkamyl1n&dl=0",
    "atac_fragment": "https://www.dropbox.com/scl/fi/nexccttzlzwr3lt0oe7eq/fragments_sorted.tsv.gz?rlkey=wrajzdz0f1g5gpx4m1vwmvcfh&st=s3hbbajy&dl=0",
}
