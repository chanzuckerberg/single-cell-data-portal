import boto3
from uuid import uuid4
from collections import defaultdict
import time


def create_query(client, query_id):
    response = client.get_database(
        CatalogName="AwsDataCatalog", DatabaseName="cellxgene_portal_dataset_download_logs_db"
    )
    response = client.start_query_execution(
        QueryString="SELECT key, requestdatetime, remoteip, operation, bytessent, useragent FROM cellxgene_portal_dataset_download_logs_db.dataset_download_logs WHERE operation like 'REST.GET.OBJECT';",
        ClientRequestToken=query_id,
        QueryExecutionContext={"Database": "cellxgene_portal_dataset_download_logs_db", "Catalog": "AwsDataCatalog"},
        ResultConfiguration={
            "OutputLocation": "s3://corpora-data-prod-logs-queries",
        },
    )
    print("Sent off query to logs!")
    return response.get("QueryExecutionId")


def get_query_results(client, query_id):
    # Wait for the query results
    results_have_not_been_calculated = True
    while results_have_not_been_calculated:
        try:
            response = client.get_query_execution(QueryExecutionId=query_id)
            status = response.get("QueryExecution").get("Status").get("State")
            if status == "SUCCEEDED":
                print("Query completed executing! Now aggregating results!")
                results_have_not_been_calculated = False
        except:
            print("Wasn't able to get query information yet. Please be patient!")
            time.sleep(1)

    response = client.get_query_results(QueryExecutionId=query_id, MaxResults=123)
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
        dataset_id = row_data[0].get("VarCharValue").split("/")[0]

        # Get existing aggregated metrics, if any.
        if dataset_id not in metadata_by_dataset:
            metadata_by_dataset[dataset_id] = {
                "curl_download_count": 0,
                "browser_download_count": 0,
                "total_download_count": 0,
            }
        dataset_metrics = metadata_by_dataset[dataset_id]

        for index, metadata in enumerate(row_data):
            if index == 0:
                dataset_metrics["total_download_count"] += 1

            if index == 2:
                ip_address = ips.add(metadata.get("VarCharValue"))

            if index == 5:
                type_of_download = metadata.get("VarCharValue")
                if "curl" in type_of_download:
                    dataset_metrics["curl_download_count"] += 1
                else:
                    dataset_metrics["browser_download_count"] += 1

        metadata_by_dataset[dataset_id] = dataset_metrics

    print(f"Total number of downloads of all datasets: {total_downloads}")
    print(f"Total number of unique IP addresses: {len(ips)}")
    print()
    for dataset_id, dataset_metrics in metadata_by_dataset.items():
        print(
            f"Dataset with id {dataset_id} was downloaded a total of {dataset_metrics['total_download_count']} times. {dataset_metrics['curl_download_count']} times via curl and {dataset_metrics['browser_download_count']} times via browser."
        )


if __name__ == "__main__":
    client = boto3.client("athena", region_name="us-west-2")
    query_id = str(uuid4())

    # Create query
    query_execution_id = create_query(client, query_id)

    # Get query results
    get_query_results(client, query_execution_id)