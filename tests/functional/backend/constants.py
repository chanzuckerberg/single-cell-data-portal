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
    "https://www.dropbox.com/scl/fi/y50umqlcrbz21a6jgu99z/5_0_0_example_valid.h5ad?rlkey"
    "=s7p6ybyx082hswix26hbl11pm&dl=0"
)

VISIUM_DATASET_URI = "https://www.dropbox.com/scl/fi/lmhue0va6ihk50ivp26da/visium_small.h5ad?rlkey=n0fo4dyi1ah7ckg9kgzwlhm8s&st=59m7g97u&dl=0"
