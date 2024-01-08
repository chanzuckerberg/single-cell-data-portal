import ast
import datetime
import json
import os
import sys

import requests
import yaml
from tqdm import tqdm


def generate_scatter_plots(json1, json2):
    import matplotlib.pyplot as plt

    # Load the profiling results
    with open(json1, "r") as f1, open(json2, "r") as f2:
        metrics = json.load(f1)["profiling_metrics_keys"]
        profiling_results1 = json.load(f1)["profiles"]
        profiling_results2 = json.load(f2)["profiles"]

    # Create a new figure for each metric
    compare_options = {d["params"]["compare"] for d in profiling_results1}
    sex_criteria = [
        ast.literal_eval(i)
        for i in {str(d["params"]["payload"]["filter"]["sex_ontology_term_ids"]) for d in profiling_results1}
    ]
    num_conditions = len(compare_options) * len(sex_criteria)

    fig, axs = plt.subplots(len(metrics), num_conditions)
    fig.set_size_inches(4 * num_conditions, 12)

    # Generate scatter plots for the collected metrics
    for i, metric in enumerate(metrics):
        for j, compare_option in enumerate(compare_options):
            for k, sex in enumerate(sex_criteria):
                metric_values1 = [
                    d[metric]
                    for d in profiling_results1
                    if d["params"]["compare"] == compare_option
                    and d["params"]["payload"]["filter"]["sex_ontology_term_ids"] == sex
                ]
                metric_values2 = [
                    d[metric]
                    for d in profiling_results2
                    if d["params"]["compare"] == compare_option
                    and d["params"]["payload"]["filter"]["sex_ontology_term_ids"] == sex
                ]

                max_value = max(metric_values1 + metric_values2)
                row, col = i, j * len(sex_criteria) + k
                axs[row, col].plot([0, max_value], [0, max_value], "k:")
                axs[row, col].scatter(metric_values1, metric_values2)
                axs[row, col].set_title(f"compare: {compare_option}, sex: {sex}")
                axs[row, col].set_xlabel("Run 1")
                axs[row, col].set_ylabel("Run 2")
                axs[row, col].set_xlim(0, max_value + 0.5)
                axs[row, col].set_ylim(0, max_value + 0.5)

    # Manually add super-titles for each row
    for i, metric in enumerate(metrics):
        # Calculate y position for the super-title of each row
        y_position = 1 - (i - 0.01) / len(metrics)
        fig.text(0.5, y_position, metric, ha="center", va="center", fontsize=14, weight="bold")

    fig.tight_layout(pad=3.0)
    plt.rcParams.update({"font.size": 10})

    # Save the scatter plots to a PDF
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    fig.savefig(f"scatter_plots_{timestamp}.pdf")
    print(f"Scatter plots output to scatter_plots_{timestamp}.pdf")


# Add the root directory to the Python module search path so you can reference backend
# without needing to move this script to the root directory to run it.
scripts_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(scripts_dir)
sys.path.append(root_dir)

# The set of random genes to use for the profiling. These are genes that are present in Census.
GENE_SETS = [
    ["ENSG00000132849"],
    ["ENSG00000245552", "ENSG00000258473", "ENSG00000230084", "ENSG00000168710", "ENSG00000146147"],
    [
        "ENSG00000262429",
        "ENSG00000204055",
        "ENSG00000068912",
        "ENSG00000244055",
        "ENSG00000258345",
        "ENSG00000284695",
        "ENSG00000255231",
        "ENSG00000287101",
        "ENSG00000277007",
        "ENSG00000166478",
    ],
    [
        "ENSG00000113522",
        "ENSG00000078814",
        "ENSG00000237111",
        "ENSG00000264582",
        "ENSG00000253684",
        "ENSG00000230343",
        "ENSG00000232586",
        "ENSG00000156966",
        "ENSG00000274139",
        "ENSG00000226592",
        "ENSG00000143569",
        "ENSG00000201421",
        "ENSG00000278172",
        "ENSG00000213658",
        "ENSG00000232113",
        "ENSG00000249162",
        "ENSG00000061938",
        "ENSG00000232305",
        "ENSG00000186106",
        "ENSG00000234595",
        "ENSG00000110002",
        "ENSG00000265182",
        "ENSG00000219361",
        "ENSG00000135249",
        "ENSG00000259882",
    ],
    [
        "ENSG00000286827",
        "ENSG00000176393",
        "ENSG00000206903",
        "ENSG00000248739",
        "ENSG00000224169",
        "ENSG00000280438",
        "ENSG00000235933",
        "ENSG00000176049",
        "ENSG00000261198",
        "ENSG00000127511",
        "ENSG00000228423",
        "ENSG00000198597",
        "ENSG00000163377",
        "ENSG00000239831",
        "ENSG00000079277",
        "ENSG00000202089",
        "ENSG00000134389",
        "ENSG00000279780",
        "ENSG00000199715",
        "ENSG00000237121",
        "ENSG00000283383",
        "ENSG00000010327",
        "ENSG00000286554",
        "ENSG00000241163",
        "ENSG00000162851",
        "ENSG00000233610",
        "ENSG00000267624",
        "ENSG00000227241",
        "ENSG00000267275",
        "ENSG00000245598",
        "ENSG00000261203",
        "ENSG00000226421",
        "ENSG00000164488",
        "ENSG00000249222",
        "ENSG00000262050",
        "ENSG00000188001",
        "ENSG00000259215",
        "ENSG00000207067",
        "ENSG00000268038",
        "ENSG00000255568",
        "ENSG00000254798",
        "ENSG00000268460",
        "ENSG00000279754",
        "ENSG00000168517",
        "ENSG00000211752",
        "ENSG00000230040",
        "ENSG00000233090",
        "ENSG00000265883",
        "ENSG00000144227",
        "ENSG00000111319",
    ],
    [
        "ENSG00000180644",
        "ENSG00000256262",
        "ENSG00000112038",
        "ENSG00000279721",
        "ENSG00000222635",
        "ENSG00000278250",
        "ENSG00000171234",
        "ENSG00000120594",
        "ENSG00000249236",
        "ENSG00000134250",
        "ENSG00000270526",
        "ENSG00000272678",
        "ENSG00000287845",
        "ENSG00000258587",
        "ENSG00000102978",
        "ENSG00000232743",
        "ENSG00000262294",
        "ENSG00000253433",
        "ENSG00000149639",
        "ENSG00000232069",
        "ENSG00000230732",
        "ENSG00000279968",
        "ENSG00000101134",
        "ENSG00000199179",
        "ENSG00000284464",
        "ENSG00000287964",
        "ENSG00000241136",
        "ENSG00000170260",
        "ENSG00000287505",
        "ENSG00000271122",
        "ENSG00000269095",
        "ENSG00000213671",
        "ENSG00000145779",
        "ENSG00000197860",
        "ENSG00000132507",
        "ENSG00000101890",
        "ENSG00000224137",
        "ENSG00000225864",
        "ENSG00000215973",
        "ENSG00000250740",
        "ENSG00000254858",
        "ENSG00000183682",
        "ENSG00000231852",
        "ENSG00000274467",
        "ENSG00000107020",
        "ENSG00000264300",
        "ENSG00000223383",
        "ENSG00000228670",
        "ENSG00000163428",
        "ENSG00000229345",
        "ENSG00000225963",
        "ENSG00000262000",
        "ENSG00000239211",
        "ENSG00000259078",
        "ENSG00000232381",
        "ENSG00000236590",
        "ENSG00000232132",
        "ENSG00000174500",
        "ENSG00000267302",
        "ENSG00000229032",
        "ENSG00000222432",
        "ENSG00000215455",
        "ENSG00000274529",
        "ENSG00000240868",
        "ENSG00000226516",
        "ENSG00000217680",
        "ENSG00000136021",
        "ENSG00000257142",
        "ENSG00000249792",
        "ENSG00000277422",
        "ENSG00000263020",
        "ENSG00000259977",
        "ENSG00000232470",
        "ENSG00000239595",
        "ENSG00000186350",
        "ENSG00000257342",
        "ENSG00000258757",
        "ENSG00000198369",
        "ENSG00000250438",
        "ENSG00000197157",
        "ENSG00000286424",
        "ENSG00000279368",
        "ENSG00000267834",
        "ENSG00000258627",
        "ENSG00000279357",
        "ENSG00000225923",
        "ENSG00000224440",
        "ENSG00000255860",
        "ENSG00000109220",
        "ENSG00000170827",
        "ENSG00000267661",
        "ENSG00000259525",
        "ENSG00000225787",
        "ENSG00000199695",
        "ENSG00000127666",
        "ENSG00000176695",
        "ENSG00000223811",
        "ENSG00000187870",
        "ENSG00000207215",
        "ENSG00000277640",
    ],
    [
        "ENSG00000277438",
        "ENSG00000225625",
        "ENSG00000222285",
        "ENSG00000276949",
        "ENSG00000244675",
        "ENSG00000207405",
        "ENSG00000147996",
        "ENSG00000204393",
        "ENSG00000120498",
        "ENSG00000234458",
        "ENSG00000270324",
        "ENSG00000163463",
        "ENSG00000221514",
        "ENSG00000213640",
        "ENSG00000284116",
        "ENSG00000218261",
        "ENSG00000231728",
        "ENSG00000202406",
        "ENSG00000255367",
        "ENSG00000176907",
        "ENSG00000248643",
        "ENSG00000234262",
        "ENSG00000112167",
        "ENSG00000224042",
        "ENSG00000187144",
        "ENSG00000102057",
        "ENSG00000200422",
        "ENSG00000262519",
        "ENSG00000207737",
        "ENSG00000259935",
        "ENSG00000279670",
        "ENSG00000125910",
        "ENSG00000236775",
        "ENSG00000213368",
        "ENSG00000234620",
        "ENSG00000286310",
        "ENSG00000237446",
        "ENSG00000216480",
        "ENSG00000169758",
        "ENSG00000286731",
        "ENSG00000258360",
        "ENSG00000183747",
        "ENSG00000229051",
        "ENSG00000273139",
        "ENSG00000272343",
        "ENSG00000200138",
        "ENSG00000182220",
        "ENSG00000224719",
        "ENSG00000091140",
        "ENSG00000105255",
        "ENSG00000111670",
        "ENSG00000262529",
        "ENSG00000207952",
        "ENSG00000287393",
        "ENSG00000185721",
        "ENSG00000237410",
        "ENSG00000238032",
        "ENSG00000109332",
        "ENSG00000225992",
        "ENSG00000183248",
        "ENSG00000233118",
        "ENSG00000279187",
        "ENSG00000274256",
        "ENSG00000119596",
        "ENSG00000223215",
        "ENSG00000200795",
        "ENSG00000282855",
        "ENSG00000081803",
        "ENSG00000143367",
        "ENSG00000252963",
        "ENSG00000233955",
        "ENSG00000214534",
        "ENSG00000275278",
        "ENSG00000252530",
        "ENSG00000257512",
        "ENSG00000255561",
        "ENSG00000254315",
        "ENSG00000264217",
        "ENSG00000224482",
        "ENSG00000238793",
        "ENSG00000263512",
        "ENSG00000166825",
        "ENSG00000260452",
        "ENSG00000201140",
        "ENSG00000259358",
        "ENSG00000186417",
        "ENSG00000143252",
        "ENSG00000223637",
        "ENSG00000202051",
        "ENSG00000197776",
        "ENSG00000072163",
        "ENSG00000225179",
        "ENSG00000224614",
        "ENSG00000231305",
        "ENSG00000200021",
        "ENSG00000233724",
        "ENSG00000264078",
        "ENSG00000131188",
        "ENSG00000258741",
        "ENSG00000234685",
        "ENSG00000204278",
        "ENSG00000143801",
        "ENSG00000251205",
        "ENSG00000261127",
        "ENSG00000105738",
        "ENSG00000227615",
        "ENSG00000207281",
        "ENSG00000155052",
        "ENSG00000286004",
        "ENSG00000277818",
        "ENSG00000260274",
        "ENSG00000253070",
        "ENSG00000239797",
        "ENSG00000137193",
        "ENSG00000238152",
        "ENSG00000138669",
        "ENSG00000209480",
        "ENSG00000233623",
        "ENSG00000243150",
        "ENSG00000151012",
        "ENSG00000249635",
        "ENSG00000265660",
        "ENSG00000185897",
        "ENSG00000197417",
        "ENSG00000237282",
        "ENSG00000228046",
        "ENSG00000226986",
        "ENSG00000286848",
        "ENSG00000250769",
        "ENSG00000153015",
        "ENSG00000224607",
        "ENSG00000101935",
        "ENSG00000178596",
        "ENSG00000267390",
        "ENSG00000287859",
        "ENSG00000176971",
        "ENSG00000228292",
        "ENSG00000213424",
        "ENSG00000256752",
        "ENSG00000249441",
        "ENSG00000275592",
        "ENSG00000156096",
        "ENSG00000145198",
        "ENSG00000104140",
        "ENSG00000213964",
        "ENSG00000101082",
        "ENSG00000255047",
        "ENSG00000172671",
        "ENSG00000211663",
        "ENSG00000244402",
        "ENSG00000240053",
        "ENSG00000237090",
        "ENSG00000214253",
        "ENSG00000275572",
        "ENSG00000156076",
        "ENSG00000136045",
        "ENSG00000282458",
        "ENSG00000279094",
        "ENSG00000287350",
        "ENSG00000287310",
        "ENSG00000233324",
        "ENSG00000229853",
        "ENSG00000212658",
        "ENSG00000183386",
        "ENSG00000255308",
        "ENSG00000106336",
        "ENSG00000228701",
        "ENSG00000261368",
        "ENSG00000136842",
        "ENSG00000252776",
        "ENSG00000233347",
        "ENSG00000167978",
        "ENSG00000201050",
        "ENSG00000261243",
        "ENSG00000224257",
        "ENSG00000277249",
        "ENSG00000275771",
        "ENSG00000251219",
        "ENSG00000136848",
        "ENSG00000277869",
        "ENSG00000146067",
        "ENSG00000204237",
        "ENSG00000238917",
        "ENSG00000115841",
        "ENSG00000266088",
        "ENSG00000283505",
        "ENSG00000198951",
        "ENSG00000146540",
        "ENSG00000162613",
        "ENSG00000267225",
        "ENSG00000231793",
        "ENSG00000230519",
        "ENSG00000105352",
        "ENSG00000260545",
        "ENSG00000266282",
        "ENSG00000260680",
        "ENSG00000238344",
        "ENSG00000227261",
        "ENSG00000254358",
        "ENSG00000237711",
    ],
]
# The groupby options to use for each query
GROUPBY_OPTIONS = [None, "disease", "publication"]
SEX_CRITERIA = [[], ["PATO:0000383"]]  # PATO:0000383 is the term id for female

# Construct the POST bodies for each combination of groupby and number of genes
POST_BODIES = []
for groupby in GROUPBY_OPTIONS:
    for sex_criteria in SEX_CRITERIA:
        for gene_set in GENE_SETS:
            payload = {
                "filter": {
                    "dataset_ids": [],
                    "development_stage_ontology_term_ids": [],
                    "disease_ontology_term_ids": [],
                    "gene_ontology_term_ids": gene_set,
                    "organism_ontology_term_id": "NCBITaxon:9606",
                    "self_reported_ethnicity_ontology_term_ids": [],
                    "sex_ontology_term_ids": sex_criteria,
                    "publication_citations": [],
                },
                "is_rollup": True,
            }
            if groupby is not None:
                payload["compare"] = groupby

            POST_BODIES.append(payload)


if __name__ == "__main__":
    with open(f"{scripts_dir}/performance_profiling_config.yaml", "r") as file:
        data = yaml.safe_load(file)

    # Parsing parameters from the config file
    rdev_auth_cookie = data["params"]["RDEV_AUTH_COOKIE"]
    api_url = data["params"]["API_URL"]
    api_url_reference = data["params"]["API_URL_REFERENCE"]
    verbose = data["params"]["VERBOSE"]
    generate_plots = data["params"]["PLOTS"]

    api_urls = [api_url]
    if api_url_reference != "":
        api_urls.append(api_url_reference)

    output_files = []
    for api_url in api_urls:
        print("Profiling API:", api_url)

        profiling_dicts = []
        profiling_metrics_keys = []
        for i, body in enumerate(tqdm(POST_BODIES, disable=verbose)):
            if verbose:
                print(f"Executing query {i+1}/{len(POST_BODIES)} with body: {body}")

            headers = {"Content-Type": "application/json"}
            if rdev_auth_cookie != "_oauth2_proxy=..." and rdev_auth_cookie != "":
                headers["Cookie"] = rdev_auth_cookie
            response = requests.post(f"{api_url}/wmg/v2/query", data=json.dumps(body), headers=headers)

            server_timing_header = response.headers.get("Server-Timing")
            assert response.status_code == 200
            assert server_timing_header is not None

            # Parse the server timing header to get the profiling results
            server_timing_values = server_timing_header.split(", ")
            profiling_dict = {}
            total_time = 0
            for value in server_timing_values:
                name, dur = value.split(";dur=")
                dur, _ = dur.split(";desc")
                dur = float(dur) / 1000
                profiling_dict[name] = dur
                total_time += dur
                profiling_metrics_keys.append(name)
            response_time = response.elapsed.total_seconds()
            download_time = response_time - total_time
            profiling_dict["total_time_backend"] = total_time
            profiling_dict["total_time_response"] = response_time
            profiling_dict["download_time"] = download_time
            response_size = len(response.content) / (1024 * 1024)
            profiling_dict["response_size_mb"] = response_size

            if verbose:
                print("Profiling results:", profiling_dict)

            profiling_dict["params"] = {
                "payload": body,
                "num_genes": len(body["filter"]["gene_ontology_term_ids"]),
                "compare": body.get("compare"),
            }
            profiling_dicts.append(profiling_dict)

        profiling_results = {
            "profiles": profiling_dicts,
            "api_url": api_url,
            "total_time_backend": sum([d["total_time_backend"] for d in profiling_dicts]),
            "total_time_response": sum([d["total_time_response"] for d in profiling_dicts]),
            "total_time_download": sum([d["download_time"] for d in profiling_dicts]),
            "total_response_size_mb": sum([d["response_size_mb"] for d in profiling_dicts]),
            "profiling_metrics_keys": profiling_metrics_keys,
        }
        print(f"Total number of queries executed: {len(profiling_dicts)}")
        print("Total time backend:", profiling_results["total_time_backend"])
        print("Total time response:", profiling_results["total_time_response"])
        print("Total time download:", profiling_results["total_time_download"])
        print("Total response size (MB):", profiling_results["total_response_size_mb"])

        # Write the profiling results to a file
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        with open(f"profiling_results_{timestamp}.json", "w") as f:
            json.dump(profiling_results, f)

            print(f"Profiling results output to {f.name}")
            output_files.append(f.name)

        # Create a new figure for each metric
        if generate_plots:
            import matplotlib.pyplot as plt

            # Generate matplotlib plots for the collected metrics
            metrics = profiling_metrics_keys + ["total_time_backend"]
            num_genes = [d["params"]["num_genes"] for d in profiling_dicts]
            compare_options = {d["params"]["compare"] for d in profiling_dicts}

            fig, axs = plt.subplots(len(metrics), 1)
            fig.set_size_inches(8, 12)
            for i, metric in enumerate(metrics):
                for compare_option in compare_options:
                    for sex_criteria in SEX_CRITERIA:
                        compare_dicts = [
                            d
                            for d in profiling_dicts
                            if d["params"]["compare"] == compare_option
                            and d["params"]["payload"]["filter"]["sex_ontology_term_ids"] == sex_criteria
                        ]
                        metric_values = [d[metric] for d in compare_dicts]
                        num_genes_compare = [d["params"]["num_genes"] for d in compare_dicts]

                        axs[i].plot(num_genes_compare, metric_values, label=f"{compare_option} {sex_criteria}")
                        axs[i].set_title(f"{metric} vs number of genes")
                        axs[i].set_xlabel("number of genes")
                        axs[i].set_ylabel(metric)
                axs[i].legend(loc="center left", bbox_to_anchor=(1, 0.5))
            fig.subplots_adjust(right=0.85)  # Adjust the whitespace on the right
            fig.tight_layout(pad=3.0)
            plt.rcParams.update({"font.size": 10})

            # Save the plots to a PDF
            fig.savefig(f"profiling_plots_{timestamp}.pdf")
            print(f"Profiling plots output to profiling_plots_{timestamp}.pdf")

    if len(api_urls) == 2:
        generate_scatter_plots(output_files[0], output_files[1])
