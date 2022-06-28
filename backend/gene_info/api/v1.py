import string
from flask import jsonify, make_response
from backend.gene_info.api.ncbi_provider import NCBIAPIException, NCBIProvider, NCBIUnexpectedResultException


def gene_info(geneID: string):
    provider = NCBIProvider()

    # search for gene UID from ensembl ID
    try:
        uid = provider.fetch_gene_uid(geneID)
    except NCBIAPIException:
        return make_response(jsonify("Failed search of NCBI database, API key issue"), 404)
    except NCBIUnexpectedResultException:
        return make_response(jsonify("Unexpected NCBI search result"), 404)

    # fetch gene information using NCBI UID
    try:
        gene_info_tree = provider.fetch_gene_info_tree(uid)
    except NCBIAPIException:
        return make_response(jsonify("Failed fetch of NCBI database, API key issue"), 404)
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
