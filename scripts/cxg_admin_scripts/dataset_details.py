import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.common.entities import Collection, Dataset, DatasetAsset
from backend.common.utils.db_session import db_session_manager


def get_public_dataset_details():
    # id, name, organism, tissue, assay, sex, cell_count, explorer_url, and S3 uris
    with db_session_manager() as session:

        published_datasets = (
            session.query(
                Dataset.table.id,
                Dataset.table.name,
                Dataset.table.organism,
                Dataset.table.tissue,
                Dataset.table.assay,
                Dataset.table.sex,
                Dataset.table.cell_count,
                Dataset.table.explorer_url,
            )
            .join(Dataset.table.collection)
            .filter(
                Dataset.table.published == "TRUE",
                Collection.table.visibility == "PUBLIC",
                Dataset.table.tombstone == "FALSE",
            )
            .all()
        )
        datasets = {}
        for dataset in published_datasets:
            dataset_data = {}
            dataset_data["name"] = dataset[1]
            dataset_data["organisms"] = [x["label"] for x in dataset[2]]
            dataset_data["tissue"] = [x["label"] for x in dataset[3]]
            dataset_data["assay"] = [x["label"] for x in dataset[4]]
            dataset_data["sex"] = [x["label"] for x in dataset[5]]
            dataset_data["cell_count"] = dataset[6]
            dataset_data["explorer_url"] = dataset[7]
            dataset_data["s3_uris"] = []
            datasets[dataset[0]] = dataset_data

        s3_uris = session.query(DatasetAsset.table.dataset_id, DatasetAsset.table.s3_uri).filter(
            DatasetAsset.table.dataset_id.in_(datasets.keys())
        )

        for s3_uri in s3_uris:
            datasets[s3_uri[0]]["s3_uris"].append(s3_uri[1])

    return datasets
