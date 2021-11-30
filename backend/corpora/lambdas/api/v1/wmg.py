import boto3
from flask import make_response


def __request_wmg_s3_object(key: str) -> str:
    s3 = boto3.client("s3")
    return (
        s3.get_object(Bucket="wmg-prototype-data-dev-public", Key=f"lung-tissue-10x-human/{key}")["Body"]
        .read()
        .decode("utf-8")
    )


def get_genes():
    return make_response(__request_wmg_s3_object("lung_tissue_genes.json"), 200)


def get_cell_types():
    return make_response(__request_wmg_s3_object("lung_tissue_cell_types.json"), 200)


def get_gene_expression(gene_name):
    return make_response(__request_wmg_s3_object(f"genes/{gene_name}.json"), 200)
