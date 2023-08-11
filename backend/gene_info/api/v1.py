import string

from flask import jsonify, make_response

from backend.common.utils.http_exceptions import ForbiddenHTTPException, NotFoundHTTPException
from backend.gene_info.api.ensembl_ids import GeneChecker
from backend.gene_info.api.ncbi_provider import NCBIAPIException, NCBIProvider, NCBIUnexpectedResultException


def gene_info(gene: string = None, geneID: string = None):  # type: ignore
    provider = NCBIProvider()

    # given just a gene name (finds corresponding gene ID)
    # in the event of a mismatch between name and ID, ID is initially selected
    # however, note that if ID fails to return a result, name will be used for ncbi search
    if not geneID:
        gene_checker = GeneChecker()
        try:
            geneID = gene_checker.get_id(gene)
        except ValueError as e:
            raise NotFoundHTTPException(str(e)) from None

    # search for gene UID from ensembl ID
    try:
        (uid, show_warning_banner) = provider.fetch_gene_uid(geneID, gene)
    except NCBIAPIException:
        raise ForbiddenHTTPException("Failed search of NCBI database, API key issue") from None
    except NCBIUnexpectedResultException:
        raise NotFoundHTTPException("Unexpected NCBI search result") from None

    # fetch gene information using NCBI UID
    try:
        gene_info_tree = provider.fetch_gene_info_tree(uid)
    except NCBIAPIException:
        raise ForbiddenHTTPException("Failed fetch of NCBI database, API key issue") from None
    except NCBIUnexpectedResultException:
        raise NotFoundHTTPException("Unexpected NCBI info tree result") from None
    gene_info = provider.parse_gene_info_tree(gene_info_tree)
    return make_response(
        jsonify(
            dict(
                name=gene_info["name"],
                summary=gene_info["summary"],
                ncbi_url=f"https://www.ncbi.nlm.nih.gov/gene/{uid}",
                synonyms=gene_info["synonyms"],
                show_warning_banner=show_warning_banner,
            )
        ),
        200,
    )
