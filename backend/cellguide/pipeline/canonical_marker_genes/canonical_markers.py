asctb_data = requests.get("https://purl.org/ccf/releases/2.2.1/ccf-asctb-all.json").json()

wmg_tissues = [
    i.split(" ")[0]
    for i in list(snapshot.cell_counts_cube.df[:]["tissue_ontology_term_id"].unique())
    if i.startswith("UBERON:")
]


def get_gene_name(gene):
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
        except:
            try:
                title = data["message"]["items"][0]["title"][0]
                citation = format_citation_mg(data["message"]["items"][0])
            except:
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


def get_tissue_name(t):
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


# get formatted citations per doi
doi_to_citation = {}
for tissue in asctb_data:
    print(tissue)
    data = asctb_data[tissue]["data"]
    for row in data:
        cell_types = [i["id"] for i in row["cell_types"] if i["id"].startswith("CL:")]
        if len(cell_types) and len(row["biomarkers_gene"]):
            refs = []
            titles = []
            for ref in row["references"]:
                doi = ref["doi"]
                if doi != "" and doi[-1] == ".":
                    doi = doi[:-1]
                if " " in doi:
                    doi = doi.split(" ")[1]
                doi = doi.strip()

                if doi != "" and doi not in doi_to_citation:
                    title = get_title_and_citation_from_doi(doi)
                    doi_to_citation[doi] = title


wmg_tissues = [
    i.split(" ")[0]
    for i in list(snapshot.cell_counts_cube.df[:]["tissue_ontology_term_id"].unique())
    if i.startswith("UBERON:")
]
all_human_genes = [list(i.values())[0] for i in snapshot.primary_filter_dimensions["gene_terms"]["NCBITaxon:9606"]]

gene_name_memory = {}

seen = []
parsed_table_entries = []

for tissue in asctb_data:
    print(tissue)
    data = asctb_data[tissue]["data"]
    for row in data:
        cell_types = [i["id"] for i in row["cell_types"] if i["id"].startswith("CL:")]
        if len(cell_types) and len(row["biomarkers_gene"]):

            # get tissue id
            tissue_id = row["anatomical_structures"][0]["id"]
            for entry in row["anatomical_structures"]:
                if entry["id"] in wmg_tissues:
                    tissue_id = entry["id"]

            # get gene symbols and names
            gene_symbols = []
            gene_names = []
            for gene in row["biomarkers_gene"]:
                symbol = gene["name"].upper()
                if symbol != "" and symbol in all_human_genes:
                    name = gene["rdfs_label"]
                    if name == "":
                        if symbol in gene_name_memory:
                            name = gene_name_memory[symbol]
                        else:
                            name = get_gene_name(symbol)
                            gene_name_memory[symbol] = name

                    gene_symbols.append(symbol)
                    gene_names.append(name)

            # get references
            refs = []
            titles = []
            for ref in row["references"]:
                doi = ref["doi"]
                if doi != "" and doi[-1] == ".":
                    doi = doi[:-1]
                if " " in doi:
                    doi = doi.split(" ")[1]
                doi = doi.strip()
                if doi != "":
                    title = doi_to_citation[doi]
                    refs.append(doi)
                    titles.append(title)

            refs = ";;".join(refs)
            titles = ";;".join(titles)
            for cell_type in cell_types:
                for i in range(len(gene_symbols)):
                    symbol = gene_symbols[i]
                    name = gene_names[i]

                    entry = {
                        "tissue": tissue_id,
                        "symbol": symbol,
                        "name": name,
                        "publication": refs,
                        "publication_titles": titles,
                        "cell_type_ontology_term_id": cell_type,
                    }
                    hashed_dict = hash(json.dumps(entry))
                    if hashed_dict not in seen:
                        parsed_table_entries.append(entry)
                        seen.append(entry)

ts = list({i["tissue"] for i in parsed_table_entries})

tissues_by_id = {t: get_tissue_name(t) for t in ts}

gene_infos = {}
for entry in parsed_table_entries:
    entry = entry.copy()
    ct = entry["cell_type_ontology_term_id"]
    del entry["cell_type_ontology_term_id"]

    a = gene_infos.get(ct, [])
    entry["tissue"] = tissues_by_id.get(entry["tissue"], entry["tissue"])
    a.append(entry)
    gene_infos[ct] = a


def func1(x):
    y = [y for y in x.values if y != ""]
    z = []
    for i in y:
        if i not in z:
            z.append(i)
    res = ";;".join(z)
    return res


def func2(x):
    res = x.values
    # pick a name that is not a gene symbol if possible
    res2 = [i not in gi2["symbol"].values for i in res]
    index = 0
    try:
        index = res2.index(True)
    except:
        pass
    return res[index]


for key in gene_infos:
    gi = pd.DataFrame(gene_infos[key])
    gi2 = pd.DataFrame(gi)

    gi = (
        gi.groupby(["tissue", "symbol"])
        .agg({"name": func2, "publication": func1, "publication_titles": func1})
        .reset_index()
        .to_dict(orient="records")
    )
    gi2 = pd.DataFrame(gi)
    #     gi2['n']=1
    #     gi3 = gi2.groupby(['symbol','tissue']).sum(numeric_only=True)['n']
    #     valid_genes = list(set(gi3.index[gi3>0].get_level_values('symbol')))
    #     gi2 = pd.DataFrame(gi)
    #     gi2 = gi2[gi2['symbol'].isin(valid_genes)]

    gi3 = gi2.groupby("symbol").agg({"name": func2, "publication": func1, "publication_titles": func1})
    gi3 = gi3.reset_index()
    gi3["tissue"] = "All Tissues"
    gi3 = gi3[gi2.columns]
    gi.extend(gi3.to_dict(orient="records"))
    gene_infos[key] = gi

json.dump(gene_infos, open("frontend/src/views/CellGuide/common/fixtures/allCellTypeMarkerGenes.json", "w"))
