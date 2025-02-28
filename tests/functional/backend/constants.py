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
    "https://www.dropbox.com/scl/fi/07l2gze8mkqy3vgh0ffh5/visium_small_5_3.h5ad?rlkey"
    "=hnktgfpce9qwpx1erqob8as2f&st=bmas6s9c&dl=0"
)

ATAC_SEQ_MANIFEST = {
    "anndata": "https://www.dropbox.com/scl/fi/erf0eglkh1pvkt4kf57xk/small_atac_seq.h5ad?rlkey=9g6sbcyl5gv83189i2i349n2p&dl=0",
    "atac_fragment": "https://www.dropbox.com/scl/fi/p4kmriyki1xyvcc9bvwxc/fragments_sorted.tsv.gz?rlkey=hydxliidfy4yneaan2rrw2arp&dl=0",
}
