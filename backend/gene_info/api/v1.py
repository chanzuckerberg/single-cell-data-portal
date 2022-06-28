import string
from flask import jsonify, make_response
from backend.gene_info.api.ncbi_provider import NCBIAPIException, NCBIProvider, NCBIUnexpectedResultException
from ...corpora.common.utils.http_exceptions import (
    NotFoundHTTPException,
    ForbiddenHTTPException,
)
import logging


def gene_info(geneID: string):
    provider = NCBIProvider()

    # search for gene UID from ensembl ID
    try:
        uid = provider.fetch_gene_uid(geneID)
    except NCBIAPIException:
        logging.error("Failed search of NCBI database, API key issue")
        raise ForbiddenHTTPException
    except NCBIUnexpectedResultException:
        logging.error("Unexpected NCBI search result")
        raise NotFoundHTTPException

    # fetch gene information using NCBI UID
    try:
        gene_info_tree = provider.fetch_gene_info_tree(uid)
    except NCBIAPIException:
        logging.error("Failed fetch of NCBI database, API key issue")
        raise ForbiddenHTTPException
    gene_info = provider.parse_gene_info_tree(gene_info_tree)
    return make_response(
        jsonify(
            dict(
                name=gene_info["name"],
                summary=gene_info["summary"],
                ncbi_url=f"https://www.ncbi.nlm.nih.gov/gene/{uid}",
                synonyms=gene_info["synonyms"],
            )
        ),
        200,
    )
