import string
from flask import jsonify, make_response
from backend.gene_info.api.ncbi_provider import NCBIAPIException, NCBIProvider, NCBIUnexpectedResultException
from backend.corpora.common.utils.http_exceptions import (
    NotFoundHTTPException,
    ForbiddenHTTPException,
)


def gene_info(gene: string, geneID: string):
    provider = NCBIProvider()

    # search for gene UID from ensembl ID
    try:
        uid = provider.fetch_gene_uid(gene, geneID)
    except NCBIAPIException:
        raise ForbiddenHTTPException("Failed search of NCBI database, API key issue")
    except NCBIUnexpectedResultException:
        raise NotFoundHTTPException("Unexpected NCBI search result")

    # fetch gene information using NCBI UID
    try:
        gene_info_tree = provider.fetch_gene_info_tree(uid)
    except NCBIAPIException:
        raise ForbiddenHTTPException("Failed fetch of NCBI database, API key issue")
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
