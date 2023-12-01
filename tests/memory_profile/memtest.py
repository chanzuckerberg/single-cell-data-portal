"""
This Python script orchestrates the download, memory metric computation, and aggregation of memory statistics for
multiple datasets. It utilizes parallel processing to efficiently handle dataset downloads, calculates memory metrics
for H5AD files using AnnData, and stores the results in a consolidated 'memory_metrics.csv' file. The script checks
existing metrics to skip processing already evaluated datasets.
"""

import os
from typing import Dict, Optional, Tuple

import anndata
import numpy as np
import pandas as pd
import requests
from scipy import sparse

POOL_SIZE = 10  # number of parallel processes to use
GB = 1024**3  # 1 Gibibyte

"""
The following DATASET variable was generated using the notebook
https://github.com/chanzuckerberg/single-cell-curation/blob/main/notebooks/curation_api/python/get_datasets.ipynb.
The output was filtered using the following code:

>>> GB = 1024**3
>>> DATASETS = []
>>> for d in datasets:
        h5ad = [a for a in d["assets"] if a["filetype"]=='H5AD'][0]
        h5ad["filesize"] /= GB
        DATASETS.append(dict(**h5ad, dataset_version_id=d["dataset_version_id"], dataset_id= d["dataset_id"]))
>>> DATASETS.sort(key=lambda x: x["filesize"])
>>> print(DATASETS)
"""
DATASETS = [
    {
        "filesize": 0.003608851693570614,
        "filetype": "H5AD",
        "url": "https://datasets.cellxgene.cziscience.com/c12960b6-465f-4769-9d61-d28a60f3080e.h5ad",
        "dataset_version_id": "c12960b6-465f-4769-9d61-d28a60f3080e",
        "dataset_id": "add5eb84-5fc9-4f01-982e-a346dd42ee82",
    },
    {
        "filesize": 0.003636792302131653,
        "filetype": "H5AD",
        "url": "https://datasets.cellxgene.cziscience.com/dae3f323-69b6-407d-9336-08a0ec0103fa.h5ad",
        "dataset_version_id": "dae3f323-69b6-407d-9336-08a0ec0103fa",
        "dataset_id": "4eb29386-de81-452f-b3c0-e00844e8c7fd",
    },
    {
        "filesize": 0.003687310963869095,
        "filetype": "H5AD",
        "url": "https://datasets.cellxgene.cziscience.com/62d2bb9a-4cc3-4765-8aa1-6bdb9216fb08.h5ad",
        "dataset_version_id": "62d2bb9a-4cc3-4765-8aa1-6bdb9216fb08",
        "dataset_id": "78d59e4a-82eb-4a61-a1dc-da974d7ea54b",
    },
]

# download a file using requests
def download_file(url: str, local_path: str) -> Tuple[str, Optional[requests.Response]]:
    """
    Download a file from a URL to a local path.
    :param url: the URL to download from
    :param local_path: the local path to download to
    :return: the local path and the response object
    """
    response = None
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        int(response.headers.get("content-length", 0))
        block_size = 1024  # 1 Kibibyte
        with open(local_path, "wb") as fp:
            for data in response.iter_content(block_size):
                fp.write(data)
    except Exception as e:
        print(f"Error downloading file from {url}: {e}")
    finally:
        if response:
            response.close()
    return local_path, response


def memtest(file_name: str) -> Dict[str, int]:
    ad = anndata.read_h5ad(file_name)
    try:
        return {
            "size_on_disk": os.path.getsize(file_name) / GB,
            "sizeof": ad.__sizeof__() / GB,
            "estimated": ad.n_obs * ad.n_vars / GB,
            "non-zeros": count_matrix_nonzero(ad),
        }
    finally:
        ad.file.close()


def count_matrix_nonzero(adata: anndata.AnnData) -> int:
    matrix = adata.X
    if adata.n_obs == 0 or adata.n_vars == 0:
        matrix_format = "dense"
    else:
        matrix_slice = matrix[0:1, 0:1]
        if isinstance(matrix_slice, sparse.spmatrix):
            matrix_format = matrix_slice.format
        elif isinstance(matrix_slice, np.ndarray):
            matrix_format = "dense"

    def chunk_matrix(
        chunk_size: int = 10_000,
    ):
        start = 0
        n = matrix.shape[0]
        for _i in range(int(n // chunk_size)):
            end = start + chunk_size
            yield matrix[start:end]
            start = end
        if start < n:
            yield matrix[start:n]

    nnz = 0
    for matrix_chunk in chunk_matrix():
        nnz += matrix_chunk.count_nonzero() if matrix_format != "dense" else np.count_nonzero(matrix_chunk)
    return nnz


def update_dataframe(dataset):
    print(f"Updating dataframe with {dataset['dataset_id']}")
    results.loc[len(results.index)] = [
        dataset["dataset_id"],
        dataset["dataset_version_id"],
        dataset["size_on_disk"],
        dataset["sizeof"],
        dataset["estimated"],
        dataset["non-zeros"],
    ]
    # write the dataframe to a csv using the PID to make it unique
    results.to_csv(f"memory_metrics_{os.getpid()}.csv", index=False)


def job(dataset):
    file_name = f"{file_path}/{dataset['dataset_id']}.h5ad"
    download_file(dataset["url"], file_name)
    dataset.update(**memtest(file_name))
    update_dataframe(dataset)
    os.remove(file_name)


# remove and pandas column


def combine_csvs():
    df = make_df()
    for file in os.listdir():
        if file.startswith("memory_metrics"):
            df1 = pd.read_csv(file)
            df = df.append(df1)
            os.remove(file)
    # use the dataset_id as the index
    df.set_index("dataset_id", inplace=True)
    df = df.drop_duplicates()
    df.to_csv("memory_metrics.csv")
    return df


def make_df():
    df = pd.DataFrame(
        columns=["dataset_id", "dataset_version_id", "size_on_disk (GB)", "sizeof (GB)", "estimated (GB)", "non-zeros"]
    )
    return df


def initialize_pool():
    global results
    results = make_df()
    global file_path
    file_path = "./datasets"


def filter_datasets():
    """skip datasets that are already in the memory_memory.csv file"""
    df = pd.read_csv("memory_metrics.csv", index_col="dataset_id")
    ds = [dataset for dataset in DATASETS if dataset["dataset_id"] not in df.index.values]
    assert ds, "No datasets to process"
    return ds


def main():
    combine_csvs()
    # datasets = filter_datasets()
    # print(f"Processing {len(datasets)} datasets")
    # with Pool(POOL_SIZE, initializer=initialize_pool) as pool:
    #     pool.map(job, datasets)
    # combine_csvs()


if __name__ == "__main__":
    main()
