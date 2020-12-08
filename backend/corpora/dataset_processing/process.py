#!/usr/bin/env python3

import os
import subprocess

import boto3
import numpy
import scanpy
import sys

from .upload import ProgressTracker, Progress, Upload
from ..common.utils import dropbox


def check_env():
    """Verify that the required environment variables are set."""
    """Verify that the required environment variables are set."""

    missing = []
    for env_var in ["DROPBOX_URL", "ARTIFACT_BUCKET", "CELLXGENE_BUCKET", "DATASET_ID"]:
        if env_var not in os.environ:
            missing.append(env_var)
    if missing:
        raise EnvironmentError(f"Missing environment variables: {missing}")


def fetch_dropbox_url(dropbox_url, local_path):
    """Given a dropbox url, download it to local_path.

    Handles fixing the url so it downloads directly.
    """

    fixed_dropbox_url = dropbox.get_download_url_from_shared_link(dropbox_url)
    if not fixed_dropbox_url:
        raise ValueError(f"Malformed Dropbox URL: {dropbox_url}")

    total_size = dropbox.get_file_info(fixed_dropbox_url)["content-length"]
    tracker = ProgressTracker(total_size)

    progress_thread = Progress(os.environ["DATASET_ID"], tracker)
    upload_thread = Upload(fixed_dropbox_url, local_path, tracker)

    progress_thread.start()
    upload_thread.start()

    upload_thread.join()
    progress_thread.join()


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

    return {
        "organism": adata.uns["organism"],
        "tissue": list(adata.obs.tissue.unique()),
        "assay": list(adata.obs.assay.unique()),
        "disease": list(adata.obs.disease.unique()),
        "sex": list(adata.obs.sex.unique()),
        "ethnicity": list(adata.obs.ethnicity.unique()),
        "development_stage": list(adata.obs.development_stage.unique()),
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

    local_filename = fetch_dropbox_url(os.environ["DROPBOX_URL"], "local.h5ad")
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

    loom_filename = make_loom(local_filename)
    cxg_dir = make_cxg(local_filename)
    seurat_filename = make_seurat(local_filename)

    s3 = boto3.client("s3")
    s3.upload_file(
        local_filename,
        os.environ["ARTIFACT_BUCKET"],
        os.path.join(os.environ["DATASET_ID"], local_filename),
        ExtraArgs={"ACL": "bucket-owner-full-control"},
    )
    s3.upload_file(
        seurat_filename,
        os.environ["ARTIFACT_BUCKET"],
        os.path.join(os.environ["DATASET_ID"], seurat_filename),
        ExtraArgs={"ACL": "bucket-owner-full-control"},
    )
    s3.upload_file(
        loom_filename,
        os.environ["ARTIFACT_BUCKET"],
        os.path.join(os.environ["DATASET_ID"], loom_filename),
        ExtraArgs={"ACL": "bucket-owner-full-control"},
    )

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


if __name__ == "__main__":
    main()
