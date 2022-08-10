# This script was created to test edge cases of the gene_info API endpoint.
# It reads from a gene ontology file and writes to a csv file successes and failures.
# See https://docs.google.com/spreadsheets/d/1ZULeY_rBsN8u4sUREmJZC3P3PL2yZm2pPrRvSCpCMHQ/edit#gid=2059332851

import csv
import urllib.request
import json

# open files
file = open('./genes_homo_sapiens_copy.csv')
counts = open('./counts.csv', 'a')
writer = csv.writer(counts)
data = csv.reader(file)

writer.writerow(['STATUS', 'ensembl ID', 'gene name', 'count', 'result'])
for row in data:
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${row[0]}&retmode=json'
    result = json.loads(urllib.request.urlopen(url).read())

    # unsuccessful search result (success is when the code returns one UID from NCBI)
    if (
        "esearchresult" not in result or 
        "idlist" not in result["esearchresult"] or 
        len(result["esearchresult"]["idlist"]) < 1
        ):
            print(f"failed on {row[0]}, {row[1]}")
            url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=${row[1]}&retmode=json'
            result = json.loads(urllib.request.urlopen(url).read())
            try:
                print(f'{len(result["esearchresult"]["idlist"])} id(s): {int(result["esearchresult"]["idlist"][0])}')
                writer.writerow(['SUCCESS ON SECOND SEARCH', row[0], row[1], len(result["esearchresult"]["idlist"]), result])
            except:
                print(f'SEARCH 2 FAILED. {result}')
                writer.writerow(['FAILED', row[0], row[1], 0, result])
    else:
        # successful search result
        print(result)
        if len(result["esearchresult"]["idlist"]) == 1:
            writer.writerow(['SUCCESS', row[0], row[1], len(result["esearchresult"]["idlist"]), result])

        # too many ID results
        else:
            writer.writerow(['TOO MANY RESULTS', row[0], row[1], len(result["esearchresult"]["idlist"]), result])

file.close()
counts.close()