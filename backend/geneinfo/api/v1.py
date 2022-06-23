import string
from flask import jsonify, make_response
import urllib.request
import json
import xml.etree.ElementTree as ET
from backend.geneinfo.config import GeneInfoConfig


def geneinfo(geneID: string):
    # load api key from secrets
    gene_info_config = GeneInfoConfig()
    gene_info_config.load()
    api_key = gene_info_config.ncbi_api_key

    # search for gene UID from ensembl ID
    search_url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                  f"db=gene&term={geneID}&api_key={api_key}&retmode=json")
    try:
        search_data = urllib.request.urlopen(search_url).read()
    except Exception:
        return make_response(jsonify("Failed search of NCBI database, API key issue"), 404)
    search_result = json.loads(search_data)
    if ("esearchresult" not in search_result or
            "idlist" not in search_result["esearchresult"] or
            len(search_result["esearchresult"]["idlist"]) < 1 or
            int(search_result["esearchresult"]["count"]) != 1):
        return make_response(jsonify("Unexpected NCBI search result"), 404)
    uid = int(search_result["esearchresult"]["idlist"][0])

    # fetch gene information using NCBI UID
    fetch_url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
                 f"db=gene&id={uid}&api_key={api_key}&retmode=xml")
    try:
        fetch_response = urllib.request.urlopen(fetch_url).read()
    except Exception:
        return make_response(jsonify("Failed fetch of NCBI gene information"), 404)

    # parse tree result
    result_tree = ET.ElementTree(ET.fromstring(fetch_response))
    root = result_tree.getroot()
    synonyms = []
    summary = ""
    name = ""
    for x in root[0]:
        if x.tag == "Entrezgene_summary":
            summary = x.text
        elif x.tag == "Entrezgene_gene":
            for y in x[0]:
                if y.tag == "Gene-ref_desc":
                    name = y.text
                elif y.tag == "Gene-ref_syn":
                    for syn in y:
                        synonyms.append(syn.text)

    return make_response(jsonify(
        dict(
            name=name,
            summary=summary,
            ncbi_url=f"https://www.ncbi.nlm.nih.gov/gene/{uid}",
            synonyms=synonyms,
        )), 200)
