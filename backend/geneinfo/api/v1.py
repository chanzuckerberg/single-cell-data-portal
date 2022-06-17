import string
from flask import jsonify
import urllib.request, json 
import xml.etree.ElementTree as ET
from backend.geneinfo.config import GeneInfoConfig

def geneinfo(geneID: string):
    gene_info_config = GeneInfoConfig()
    gene_info_config.load()
    api_key = gene_info_config.ncbi_api_key

    # search for gene UID from ensembl ID
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term={}&api_key={}&retmode=json".format(geneID, api_key)
    print("search_url: " + search_url)
    search_data = urllib.request.urlopen(search_url).read()
    search_result = json.loads(search_data)
    uid = int(search_result["esearchresult"]["idlist"][0])

    # fetch gene information using NCBI UID
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={}&api_key={}&retmode=xml".format(uid, api_key)
    fetch_response = urllib.request.urlopen(fetch_url).read()

    # parse tree result
    result_tree = ET.ElementTree(ET.fromstring(fetch_response))
    root = result_tree.getroot()
    synonyms = []
    for x in root[0]:
        if x.tag == 'Entrezgene_summary':
            summary = x.text
        elif x.tag == 'Entrezgene_gene':
            for y in x[0]:
                if y.tag == 'Gene-ref_desc':
                    name = y.text
                elif y.tag == 'Gene-ref_syn':
                    for syn in y:
                        synonyms.append(syn.text)

    return jsonify(
        dict(
            name=name,
            summary=summary,
            ncbi_url=f"https://www.ncbi.nlm.nih.gov/gene/{uid}",
            synonyms=synonyms,
        )
    )


