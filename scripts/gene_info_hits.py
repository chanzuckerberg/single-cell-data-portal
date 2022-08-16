# This script was created to test edge cases of the gene_info API endpoint.
# It reads from a gene ontology file and writes to a csv file successes and failures.
# See https://docs.google.com/spreadsheets/d/1ZULeY_rBsN8u4sUREmJZC3P3PL2yZm2pPrRvSCpCMHQ/edit#gid=2059332851

import csv
import urllib.request
import json


def get_search_url(term):
    return f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${term}&retmode=json"


def write_row(writer, status, row, result):
    writer.writerow([status, row[0], row[1], len(result["esearchresult"]["idlist"]), result])


# open files
file = open("./genes_homo_sapiens_copy.csv")
counts = open("./counts.csv", "a")
writer = csv.writer(counts)
data = csv.reader(file)

# write column names
writer.writerow(["STATUS", "ensembl ID", "gene name", "count", "result"])

for row in data:
    result = json.loads(urllib.request.urlopen(get_search_url(row[0])).read())
    # unsuccessful search result (success is when the code returns one UID from NCBI)
    if (
        "esearchresult" not in result
        or "idlist" not in result["esearchresult"]
        or len(result["esearchresult"]["idlist"]) < 1
    ):
        result = json.loads(urllib.request.urlopen(get_search_url(row[1])).read())
        try:
            write_row(writer, "SUCCESS ON SECOND SEARCH", row, result)
        except:
            write_row(writer, "FAILED", row, result)
    else:
        if len(result["esearchresult"]["idlist"]) == 1:
            write_row(writer, "SUCCESS", row, result)
        else:
            write_row(writer, "TOO MANY RESULTS", row, result)

file.close()
counts.close()
