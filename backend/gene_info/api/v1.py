import string
from flask import jsonify, make_response
import xml.etree.ElementTree as ET
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
    gene_info = parse_gene_info_tree(gene_info_tree)
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


def parse_gene_info_tree(tree_response):
    # parse NCBI XML response into relevant values to return by gene_info API
    result_tree = ET.ElementTree(ET.fromstring(tree_response))
    root = result_tree.getroot()
    synonyms = []
    summary = ""
    name = ""
    if len(root) > 0:
        for x in root[0]:
            if x.tag == "Entrezgene_summary":
                summary = x.text
            elif x.tag == "Entrezgene_gene":
                if len(x) > 0:
                    for y in x[0]:
                        if y.tag == "Gene-ref_desc":
                            name = y.text
                        elif y.tag == "Gene-ref_syn":
                            for syn in y:
                                synonyms.append(syn.text)
    return dict(
        name=name,
        summary=summary,
        synonyms=synonyms,
    )
