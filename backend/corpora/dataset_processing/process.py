#!/usr/bin/env python3

import os
import subprocess

from os.path import basename, join

import boto3
import numpy
import scanpy
import sys

<<<<<<< HEAD
<<<<<<< HEAD
try:
    from ..common.entities.dataset import Dataset
    from ..common.corpora_orm import DatasetArtifactFileType, DatasetArtifactType
    from ..common.utils import dropbox
    from ..common.utils.db_utils import db_session
    from .download import download

except ImportError:
    # We're in the container
    sys.path.append("/code")
    from common.entities.dataset import Dataset
    from common.corpora_orm import DatasetArtifactFileType, DatasetArtifactType
    from common.utils import dropbox
    from common.utils.db_utils import db_session
    from .download import download

# This is unfortunate, but this information doesn't appear to live anywhere
# accessible to the uploader
DEPLOYMENT_STAGE_TO_URL = {
    "dev": "https://cellxgene.dev.single-cell.czi.technology",
    "staging": "https://cellxgene.staging.single-cell.czi.technology",
    "prod": "https://cellxgene.cziscience.com",
}
from .upload import ProgressTracker, Progress, Upload
=======
from .upload import upload
>>>>>>> Uploading with status updates
from ..common.utils import dropbox


def check_env():
    """Verify that the required environment variables are set."""
    """Verify that the required environment variables are set."""

    missing = []
    for env_var in ["DROPBOX_URL", "ARTIFACT_BUCKET", "CELLXGENE_BUCKET", "DATASET_ID", "DEPLOYMENT_STAGE"]:
        if env_var not in os.environ:
            missing.append(env_var)
    if missing:
        raise EnvironmentError(f"Missing environment variables: {missing}")


<<<<<<< HEAD
def create_artifacts(h5ad_filename, seurat_filename, loom_filename):

    s3 = boto3.client("s3")
    artifacts = []

    s3.upload_file(
        h5ad_filename,
        os.environ["ARTIFACT_BUCKET"],
        join(os.environ["DATASET_ID"], basename(h5ad_filename)),
        ExtraArgs={"ACL": "bucket-owner-full-control"},
    )

    artifacts.append(
        {
            "filename": basename(h5ad_filename),
            "filetype": DatasetArtifactFileType.H5AD,
            "type": DatasetArtifactType.REMIX,
            "user_submitted": True,
            "s3_uri": join("s3://", os.environ["ARTIFACT_BUCKET"], os.environ["DATASET_ID"], basename(h5ad_filename)),
        }
    )

    s3.upload_file(
        seurat_filename,
        os.environ["ARTIFACT_BUCKET"],
        join(os.environ["DATASET_ID"], basename(seurat_filename)),
        ExtraArgs={"ACL": "bucket-owner-full-control"},
    )
    artifacts.append(
        {
            "filename": basename(seurat_filename),
            "filetype": DatasetArtifactFileType.RDS,
            "type": DatasetArtifactType.REMIX,
            "user_submitted": True,
            "s3_uri": join("s3://", os.environ["ARTIFACT_BUCKET"], os.environ["DATASET_ID"], basename(seurat_filename)),
        }
    )

    s3.upload_file(
        loom_filename,
        os.environ["ARTIFACT_BUCKET"],
        join(os.environ["DATASET_ID"], basename(loom_filename)),
        ExtraArgs={"ACL": "bucket-owner-full-control"},
    )
    artifacts.append(
        {
            "filename": basename(loom_filename),
            "filetype": DatasetArtifactFileType.LOOM,
            "type": DatasetArtifactType.REMIX,
            "user_submitted": True,
            "s3_uri": join("s3://", os.environ["ARTIFACT_BUCKET"], os.environ["DATASET_ID"], basename(loom_filename)),
        }
    )

    return artifacts


@db_session()
def update_db(metadata=None, processing_status=None):

    dataset = Dataset.get(os.environ["DATASET_ID"])

    if metadata:

        # TODO: Delete this line once mean_genes_per_cell is in the db
        metadata.pop("mean_genes_per_cell", None)

        dataset.update(**metadata)

    if processing_status:
        status = dataset.processing_status.to_dict()
        status.pop("dataset")
        status.pop("created_at")
        status.pop("updated_at")
        status.update(processing_status)
        dataset.update(processing_status=status)


<<<<<<< HEAD
<<<<<<< HEAD
def download_from_dropbox_url(dataset_uuid: str, dropbox_url: str, local_path: str) -> str:
=======
def fetch_dropbox_url(dataset_uuid, dropbox_url, local_path):
>>>>>>> Uploading with status updates
=======
def fetch_dropbox_url(dataset_uuid: str, dropbox_url: str, local_path: str) -> str:
>>>>>>> cleanup
=======
def download_from_dropbox_url(dataset_uuid: str, dropbox_url: str, local_path: str) -> str:
>>>>>>> Making suggest changes
    """Given a dropbox url, download it to local_path.

    Handles fixing the url so it downloads directly.
    """

    fixed_dropbox_url = dropbox.get_download_url_from_shared_link(dropbox_url)
    if not fixed_dropbox_url:
        raise ValueError(f"Malformed Dropbox URL: {dropbox_url}")

    file_info = dropbox.get_file_info(fixed_dropbox_url)
<<<<<<< HEAD
    download(dataset_uuid, fixed_dropbox_url, local_path, file_info["size"])
    return local_path
    total_size = dropbox.get_file_info(fixed_dropbox_url)["content-length"]
    tracker = ProgressTracker(total_size)

    progress_thread = Progress(os.environ["DATASET_ID"], tracker)
    upload_thread = Upload(fixed_dropbox_url, local_path, tracker)

    progress_thread.start()
    upload_thread.start()

    upload_thread.join()
    progress_thread.join()
=======
    upload(dataset_uuid, fixed_dropbox_url, local_path, file_info["size"])
<<<<<<< HEAD
>>>>>>> Uploading with status updates
=======
    return local_path
>>>>>>> cleanup


def extract_metadata(filename):
    """Pull metadata out of the AnnData file to insert into the dataset table."""

    adata = scanpy.read_h5ad(filename, backed="r")

    try:
        raw_layer_name = [k for k, v in adata.uns["layer_descriptions"].items() if v == "raw"][0]
    except (KeyError, IndexError):
        raise RuntimeError("Raw layer not found in layer descriptions!")
    if raw_layer_name == "X":
        raw_layer = adata.X
    elif raw_layer_name == "raw.X":
        raw_layer = adata.raw.X
    else:
        raw_layer = adata.layers[raw_layer]

    # Calling np.count_nonzero on and h5py.Dataset appears to read the entire thing
    # into memory, so we need to chunk it to be safe.
    stride = 50000
    numerator, denominator = 0, 0
    for bounds in zip(range(0, raw_layer.shape[0], stride), range(stride, raw_layer.shape[0] + stride, stride)):
        chunk = raw_layer[bounds[0] : bounds[1], :]
        numerator += numpy.count_nonzero(chunk)
        denominator += chunk.shape[0]

    def _get_term_pairs(base_term):
        base_term_id = base_term + "_ontology_term_id"
        return [
            {"label": k[0], "ontology_term_id": k[1]}
            for k in adata.obs.groupby([base_term, base_term_id]).groups.keys()
        ]

    return {
        "organism": {"label": adata.uns["organism"], "ontology_term_id": adata.uns["organism_ontology_term_id"]},
        "tissue": _get_term_pairs("tissue"),
        "assay": _get_term_pairs("assay"),
        "disease": _get_term_pairs("disease"),
        "sex": list(adata.obs.sex.unique()),
        "ethnicity": _get_term_pairs("ethnicity"),
        "development_stage": _get_term_pairs("development_stage"),
        "cell_count": adata.shape[0],
        "mean_genes_per_cell": numerator / denominator,
    }


def make_loom(local_filename):
    """Create a loom file from the AnnData file."""

    adata = scanpy.read_h5ad(local_filename)
    column_name_map = {}
    for column in adata.obs.columns:
        if "/" in column:
            column_name_map[column] = column.replace("/", "-")
    if column_name_map:
        adata.obs = adata.obs.rename(columns=column_name_map)

    loom_filename = local_filename.replace(".h5ad", ".loom")
    adata.write_loom(loom_filename, True)
    return loom_filename


def make_seurat(local_filename):
    """Create a Seurat rds file from the AnnData file."""

    seurat_proc = subprocess.run(["Rscript", "/dataset_processing/make_seurat.R", local_filename], capture_output=True)
    if seurat_proc.returncode != 0:
        raise RuntimeError(f"Seurat conversion failed: {seurat_proc.stdout} {seurat_proc.stderr}")

    return local_filename.replace(".h5ad", ".rds")


def make_cxg(local_filename):
    cxg_dir = local_filename.replace(".h5ad", ".cxg")
    cxg_proc = subprocess.run(
        ["cellxgene", "convert", "-o", cxg_dir, "-s", "10.0", local_filename], capture_output=True
    )
    if cxg_proc.returncode != 0:
        raise RuntimeError(f"CXG conversion failed: {cxg_proc.stderr}")
    return cxg_dir


def main():
    check_env()

<<<<<<< HEAD
<<<<<<< HEAD
    local_filename = download_from_dropbox_url(os.environ["DATASET_ID"], os.environ["DROPBOX_URL"], "local.h5ad")
=======
    local_filename = fetch_dropbox_url(os.environ["DATASET_ID"], os.environ["DROPBOX_URL"], "local.h5ad")
>>>>>>> Uploading with status updates
=======
    local_filename = download_from_dropbox_url(os.environ["DATASET_ID"], os.environ["DROPBOX_URL"], "local.h5ad")
>>>>>>> Making suggest changes
    print("Download complete", flush=True)
    val_proc = subprocess.run(["cellxgene", "schema", "validate", local_filename], capture_output=True)
    if False and val_proc.returncode != 0:
        print("Validation failed!")
        print(f"stdout: {val_proc.stdout}")
        print(f"stderr: {val_proc.stderr}")
        sys.exit(1)
    print("Validation complete", flush=True)

    metadata_dict = extract_metadata(local_filename)
    print(metadata_dict, flush=True)
    update_db(metadata=metadata_dict)

    loom_filename = make_loom(local_filename)
    cxg_dir = make_cxg(local_filename)
    seurat_filename = make_seurat(local_filename)
    artifacts = create_artifacts(local_filename, seurat_filename, loom_filename)

    subprocess.run(
        [
            "aws",
            "s3",
            "cp",
            cxg_dir,
            f"s3://{os.environ['CELLXGENE_BUCKET']}/{os.environ['DATASET_ID']}/",
            "--recursive",
            "--acl",
            "bucket-owner-full-control",
        ],
        check=True,
    )
    deployment_directories = [
        {"url": join(DEPLOYMENT_STAGE_TO_URL[os.environ["DEPLOYMENT_STAGE"]], os.environ["DATASET_ID"], "")}
    ]

    update_db(metadata={"artifacts": artifacts, "deployment_directories": deployment_directories})


if __name__ == "__main__":
    main()
