import boto3
from uuid import uuid4
from collections import defaultdict
import time
from pandas import DataFrame
import datetime
import numpy as np
import pandas
import requests

pandas.set_option(
    "display.max_rows", None, "display.max_columns", None, "display.width", 1000, "display.max_colwidth", None
)


def get_json_response(collection_id):
    collection_info_response = requests.get(f"https://api.cellxgene.cziscience.com/dp/v1/collections/{collection_id}")
    collection_info = collection_info_response.json()
    return collection_info


def construct_map_of_dataset_assets_to_collection_and_dataset_information_from_apis():
    # Fetch list of collection ids
    collections_response = requests.get("https://api.cellxgene.cziscience.com/dp/v1/collections")
    collections = collections_response.json()["collections"]
    collection_ids = [collection.get("id") for collection in collections]

    all_collection_responses = {}
    for collection_id in collection_ids:
        all_collection_responses[collection_id] = get_json_response(collection_id)

    dataset_name_by_s3_uri = {}
    collection_id_by_s3_uri = {}

    for collection_id, collection_information in all_collection_responses.items():
        datasets_in_collection = collection_information.get("datasets")

        for dataset in datasets_in_collection:
            dataset_id = dataset.get("id")
            dataset_name = dataset.get("name")
            dataset_assets = dataset.get("dataset_assets")

            for dataset_asset in dataset_assets:
                uri = dataset_asset.get("s3_uri")
                uri = uri.replace("s3://corpora-data-prod/", "")

                if uri in dataset_name_by_s3_uri:
                    print(f"ERROR: Why is this URI repeated in dataset name dict? {uri}")
                else:
                    dataset_name_by_s3_uri[uri] = dataset_name
                if uri in collection_id_by_s3_uri:
                    print(f"ERROR: Why is this URI repeated in collection id dict? {uri}")
                else:
                    collection_id_by_s3_uri[uri] = collection_id

    return dataset_name_by_s3_uri, collection_id_by_s3_uri


def create_query(client, query_id, today_datestring=None):
    if today_datestring is None:
        today = datetime.date.today().strftime("%Y-%m-%d:%H:%M:%S")
        one_week_ago = (datetime.datetime.now() - datetime.timedelta(days=7)).date().strftime("%Y-%m-%d:%H:%M:%S")
    else:
        today_datetime = datetime.date.fromisoformat(today_datestring)
        today = today_datetime.strftime("%Y-%m-%d:%H:%M:%S")
        one_week_ago = (
            (datetime.datetime.fromisoformat(today_datestring) - datetime.timedelta(days=7))
            .date()
            .strftime("%Y-%m-%d:%H:%M:%S")
        )

    print(f"Starting date is: {one_week_ago}. Ending date is: {today}.")

    query_string = (
        "SELECT key, requestdatetime, remoteip, operation, bytessent, useragent FROM "
        "cellxgene_portal_dataset_download_logs_db.dataset_download_logs WHERE operation like "
        "'REST.GET.OBJECT' AND parse_datetime(requestdatetime,'dd/MMM/yyyy:HH:mm:ss Z') BETWEEN "
        f"parse_datetime('{one_week_ago}','yyyy-MM-dd:HH:mm:ss') AND parse_datetime('{today}',"
        "'yyyy-MM-dd:HH:mm:ss');"
    )

    response = client.get_database(
        CatalogName="AwsDataCatalog", DatabaseName="cellxgene_portal_dataset_download_logs_db"
    )
    response = client.start_query_execution(
        QueryString=query_string,
        ClientRequestToken=query_id,
        QueryExecutionContext={"Database": "cellxgene_portal_dataset_download_logs_db", "Catalog": "AwsDataCatalog"},
        ResultConfiguration={
            "OutputLocation": "s3://corpora-data-prod-logs-queries",
        },
    )
    return response.get("QueryExecutionId")


def get_query_results(client, query_id, dataset_name_by_s3_uri, collection_id_by_s3_uri):
    # Wait for the query results
    results_have_not_been_calculated = True
    while results_have_not_been_calculated:
        try:
            response = client.get_query_execution(QueryExecutionId=query_id)
            status = response.get("QueryExecution").get("Status").get("State")
            if status == "SUCCEEDED":
                results_have_not_been_calculated = False
        except:
            print(f"Wasn't able to get query information for query ID {query_id} yet. Please be patient!")
            time.sleep(1)

    response = client.get_query_results(QueryExecutionId=query_id)
    rows = response.get("ResultSet").get("Rows")

    # Structure that will hold all the metrics
    metadata_by_dataset = defaultdict(dict)
    total_downloads = 0
    ips = set()

    # Delete the first row since it just contains the headers
    rows = rows[1:]

    for row in rows:
        row_data = row.get("Data")
        total_downloads += 1

        # Get dataset id
        dataset_id = row_data[0].get("VarCharValue")

        # If this dataset is private, skip its download metrics
        if (dataset_id not in collection_id_by_s3_uri) and (dataset_id not in dataset_name_by_s3_uri):
            continue

        # Generate the unique dataset index for this dataset by concatenating the collection ID and the dataset ID
        dataset_index = f"{collection_id_by_s3_uri.get(dataset_id, 'PRIVATE_COLLECTION')}/{dataset_name_by_s3_uri.get(dataset_id, 'PRIVATE_DATASET')}"

        # Get existing aggregated metrics, if any.
        if dataset_index not in metadata_by_dataset:
            metadata_by_dataset[dataset_index] = {
                "curl_downloads": 0,
                "browser_downloads": 0,
                "total_downloads": 0,
                "h5ad_downloads": 0,
                "loom_downloads": 0,
                "seurat_downloads": 0,
            }
        dataset_metrics = metadata_by_dataset[dataset_index]

        for index, metadata in enumerate(row_data):
            if index == 0:
                dataset_metrics["total_downloads"] += 1

                if "h5ad" in metadata.get("VarCharValue"):
                    dataset_metrics["h5ad_downloads"] += 1
                elif "loom" in metadata.get("VarCharValue"):
                    dataset_metrics["loom_downloads"] += 1
                elif "rds" in metadata.get("VarCharValue"):
                    dataset_metrics["seurat_downloads"] += 1

            if index == 2:
                ips.add(metadata.get("VarCharValue"))

            if index == 5:
                type_of_download = metadata.get("VarCharValue")
                if "curl" in type_of_download:
                    dataset_metrics["curl_downloads"] += 1
                else:
                    dataset_metrics["browser_downloads"] += 1

        metadata_by_dataset[dataset_index] = dataset_metrics

    dataset_metrics_df = DataFrame.from_dict(metadata_by_dataset, orient="index")
    dataset_metrics_df.to_csv("results.csv")

    print(f"Total number of downloads of all datasets: {total_downloads}")
    print(f"Total number of unique IP addresses: {len(ips)}")


if __name__ == "__main__":
    client = boto3.client("athena", region_name="us-west-2")
    query_id = str(uuid4())

    # Construct dataset id -> collection id and dataset id -> dataset name indices
    (
        dataset_name_by_s3_uri,
        collection_id_by_s3_uri,
    ) = construct_map_of_dataset_assets_to_collection_and_dataset_information_from_apis()

    # Create query
    query_execution_id = create_query(client, query_id)

    # Get query results
    get_query_results(client, query_execution_id, dataset_name_by_s3_uri, collection_id_by_s3_uri)
