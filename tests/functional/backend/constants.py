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
    "https://www.dropbox.com/scl/fi/8yizrdfcfl02dtk3ke4sg/example_5_3_valid.h5ad?rlkey"
    "=i1qc5qai9w2o9l1fithyatxdf&st=uxgudiwz&dl=0"
)

VISIUM_DATASET_URI = (
    "https://www.dropbox.com/scl/fi/3y22olsc70of8rbb1es77/visium_small.h5ad?rlkey"
    "=cgwd59ouk340zlqh6fcnthizz&st=u2nyo3xp&dl=0"
)
