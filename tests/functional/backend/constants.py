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

DATASET_URI = (
    "https://www.dropbox.com/scl/fi/593x1pa1bdnb4rwqnlp8y/example_valid.h5ad?rlkey"
    "=qzgsr0e0j5i13nsr10p983atf&st=o3b21tl4&dl=0"
)

VISIUM_DATASET_URI = (
    "https://www.dropbox.com/scl/fi/w8527i9n16461r2ydpbph/small_visium.h5ad?rlkey"
    "=x4zoxoh896o9mj2i372sxgahz&st=qijzavx6&dl=0"
)

DATASET_MANIFEST = {"anndata": DATASET_URI}
VISIUM_DATASET_MANIFEST = {"anndata": VISIUM_DATASET_URI}
ATAC_SEQ_MANIFEST = {
    "anndata": "https://www.dropbox.com/scl/fi/wxgczoo7gfu1n8fmxt350/atac.h5ad?rlkey=zbsknm1xyuzv83tixsnrbbh13&st=wbgpuzj8&dl=0",
    "atac_fragment": "https://www.dropbox.com/scl/fi/nexccttzlzwr3lt0oe7eq/fragments_sorted.tsv.gz?rlkey=wrajzdz0f1g5gpx4m1vwmvcfh&st=s3hbbajy&dl=0",
    "flags": {"deduplicate_fragments": True},
}
