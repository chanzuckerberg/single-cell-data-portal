import requests


def clean_doi(doi):
    if doi != "" and doi[-1] == ".":
        doi = doi[:-1]
    if " " in doi:
        doi = doi.split(" ")[1]
    doi = doi.strip()
    return doi


def get_gene_name_from_gene_info(gene):
    a = requests.get(f"https://api.cellxgene.dev.single-cell.czi.technology/gene_info/v1/gene_info?gene={gene}")
    if a.status_code == 200:
        r = a.json()
        return r["name"]
    else:
        return gene


def get_title_and_citation_from_doi(doi):
    url = f"https://api.crossref.org/works/{doi}"

    # Send a GET request to the API
    response = requests.get(url)

    # If the GET request is successful, the status code will be 200
    if response.status_code == 200:
        # Get the response data
        data = response.json()

        # Get the title and citation count from the data
        try:
            title = data["message"]["title"][0]
            citation = format_citation_mg(data["message"])
        except Exception:
            try:
                title = data["message"]["items"][0]["title"][0]
                citation = format_citation_mg(data["message"]["items"][0])
            except Exception:
                return doi
        return f"{title}\n\n - {citation}"
    else:
        return doi


def format_citation_mg(message):
    first_author = message["author"][0]
    if "family" in first_author:
        author_str = f"{first_author['family']}, {first_author['given']} et al."
    else:
        author_str = f"{first_author['name']} et al."

    journal = " " + message["container-title"][0] if len(message["container-title"]) else ""
    year = message["created"]["date-parts"][0][0]

    return f"{author_str} ({year}){journal}"


def get_tissue_name_from_uberon(t):
    t = t.replace(":", "_")
    urls = [
        f"https://www.ebi.ac.uk/ols4/api/ontologies/clo/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F{t}",
        f"https://www.ebi.ac.uk/ols4/api/ontologies/envo/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F{t}",
        f"https://www.ebi.ac.uk/ols4/api/ontologies/flopo/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F{t}",
        f"https://www.ebi.ac.uk/ols4/api/ontologies/doid/terms/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F{t}",
    ]
    for url in urls:
        response = requests.get(url)
        if response.status_code == 200:
            r = response.json()
            return r["label"]
    return t
