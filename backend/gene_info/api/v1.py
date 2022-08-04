import string
from flask import jsonify, make_response
from backend.gene_info.api.ncbi_provider import NCBIAPIException, NCBIProvider, NCBIUnexpectedResultException
from backend.corpora.common.utils.http_exceptions import (
    NotFoundHTTPException,
    ForbiddenHTTPException,
)
from backend.gene_info.api.ensembl_ids import GeneChecker


def gene_info(gene: string = None, geneID: string = None):
    provider = NCBIProvider()

    # given just a gene name (finds corresponding gene ID)
    # in the event of a mismatch between name and ID, ID is initially selected
    # however, note that if ID fails to return a result, name will be used for ncbi search
    if not geneID:
        gene_checker = GeneChecker()
        try:
            geneID = gene_checker.get_id(gene)
        except ValueError as e:
            raise NotFoundHTTPException(str(e))

    # search for gene UID from ensembl ID
    try:
        (uid, is_ensembl_id_result) = provider.fetch_gene_uid(geneID, gene)
    except NCBIAPIException:
        raise ForbiddenHTTPException("Failed search of NCBI database, API key issue")
    except NCBIUnexpectedResultException:
        raise NotFoundHTTPException("Unexpected NCBI search result")

    # fetch gene information using NCBI UID
    try:
        gene_info_tree = provider.fetch_gene_info_tree(uid)
    except NCBIAPIException:
        raise ForbiddenHTTPException("Failed fetch of NCBI database, API key issue")
    except NCBIUnexpectedResultException:
        raise NotFoundHTTPException("Unexpected NCBI info tree result")
    gene_info = provider.parse_gene_info_tree(gene_info_tree)
    return make_response(
        jsonify(
            dict(
                name=gene_info["name"],
                summary=gene_info["summary"],
                ncbi_url=f"https://www.ncbi.nlm.nih.gov/gene/{uid}",
                synonyms=gene_info["synonyms"],
                is_ensembl_id_result=is_ensembl_id_result,
            )
        ),
        200,
    )
