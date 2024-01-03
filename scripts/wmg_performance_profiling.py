import argparse
import datetime
import json
import os
import sys

import requests

# Add the root directory to the Python module search path so you can reference backend
# without needing to move this script to the root directory to run it.
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(root_dir)

# The set of random genes to use for the profiling. These are genes that are present in Census.
RANDOM_GENES = [
    "ENSG00000161267",
    "ENSG00000267489",
    "ENSG00000213755",
    "ENSG00000242021",
    "ENSG00000253783",
    "ENSG00000269966",
    "ENSG00000262870",
    "ENSG00000255313",
    "ENSG00000284324",
    "ENSG00000204584",
    "ENSG00000257726",
    "ENSG00000249883",
    "ENSG00000196990",
    "ENSG00000279261",
    "ENSG00000235962",
    "ENSG00000254190",
    "ENSG00000170482",
    "ENSG00000237574",
    "ENSG00000174886",
    "ENSG00000207016",
    "ENSG00000018625",
    "ENSG00000181192",
    "ENSG00000163466",
    "ENSG00000165655",
    "ENSG00000229226",
    "ENSG00000223942",
    "ENSG00000226209",
    "ENSG00000254793",
    "ENSG00000221866",
    "ENSG00000254703",
    "ENSG00000240785",
    "ENSG00000197114",
    "ENSG00000173213",
    "ENSG00000207196",
    "ENSG00000233286",
    "ENSG00000224273",
    "ENSG00000180626",
    "ENSG00000250900",
    "ENSG00000225308",
    "ENSG00000232893",
    "ENSG00000233583",
    "ENSG00000170777",
    "ENSG00000232622",
    "ENSG00000137877",
    "ENSG00000264486",
    "ENSG00000229832",
    "ENSG00000235997",
    "ENSG00000212589",
    "ENSG00000243095",
    "ENSG00000230309",
    "ENSG00000070729",
    "ENSG00000110169",
    "ENSG00000287170",
    "ENSG00000225443",
    "ENSG00000229941",
    "ENSG00000250141",
    "ENSG00000207993",
    "ENSG00000109805",
    "ENSG00000183287",
    "ENSG00000233179",
    "ENSG00000264982",
    "ENSG00000257932",
    "ENSG00000247708",
    "ENSG00000267382",
    "ENSG00000266015",
    "ENSG00000247345",
    "ENSG00000106526",
    "ENSG00000257193",
    "ENSG00000243008",
    "ENSG00000221083",
    "ENSG00000287260",
    "ENSG00000242120",
    "ENSG00000286834",
    "ENSG00000102409",
    "ENSG00000267002",
    "ENSG00000230219",
    "ENSG00000286889",
    "ENSG00000230066",
    "ENSG00000283431",
    "ENSG00000260689",
    "ENSG00000203737",
    "ENSG00000269243",
    "ENSG00000288097",
    "ENSG00000238085",
    "ENSG00000268296",
    "ENSG00000234645",
    "ENSG00000258427",
    "ENSG00000271978",
    "ENSG00000213386",
    "ENSG00000288575",
    "ENSG00000049192",
    "ENSG00000225099",
    "ENSG00000212127",
    "ENSG00000283712",
    "ENSG00000213047",
    "ENSG00000169925",
    "ENSG00000261357",
    "ENSG00000284416",
    "ENSG00000229992",
    "ENSG00000181541",
    "ENSG00000199461",
    "ENSG00000235069",
    "ENSG00000226188",
    "ENSG00000270154",
    "ENSG00000104518",
    "ENSG00000286454",
    "ENSG00000287579",
    "ENSG00000155380",
    "ENSG00000279725",
    "ENSG00000208017",
    "ENSG00000235865",
    "ENSG00000267800",
    "ENSG00000141027",
    "ENSG00000196632",
    "ENSG00000279605",
    "ENSG00000147862",
    "ENSG00000263788",
    "ENSG00000054148",
    "ENSG00000226668",
    "ENSG00000198690",
    "ENSG00000173421",
    "ENSG00000203408",
    "ENSG00000279253",
    "ENSG00000226089",
    "ENSG00000274380",
    "ENSG00000260835",
    "ENSG00000255814",
    "ENSG00000230623",
    "ENSG00000252840",
    "ENSG00000096996",
    "ENSG00000273486",
    "ENSG00000276326",
    "ENSG00000235328",
    "ENSG00000261385",
    "ENSG00000233776",
    "ENSG00000253187",
    "ENSG00000100387",
    "ENSG00000282508",
    "ENSG00000111728",
    "ENSG00000168685",
    "ENSG00000219575",
    "ENSG00000134765",
    "ENSG00000283933",
    "ENSG00000229607",
    "ENSG00000235777",
    "ENSG00000272362",
    "ENSG00000233714",
    "ENSG00000272426",
    "ENSG00000227244",
    "ENSG00000235012",
    "ENSG00000100146",
    "ENSG00000257220",
    "ENSG00000136634",
    "ENSG00000225561",
    "ENSG00000231724",
    "ENSG00000266905",
    "ENSG00000262160",
    "ENSG00000252864",
    "ENSG00000254835",
    "ENSG00000272256",
    "ENSG00000238386",
    "ENSG00000121381",
    "ENSG00000261314",
    "ENSG00000244060",
    "ENSG00000184634",
    "ENSG00000250261",
    "ENSG00000273008",
    "ENSG00000068615",
    "ENSG00000223974",
    "ENSG00000259304",
    "ENSG00000177414",
    "ENSG00000202190",
    "ENSG00000283679",
    "ENSG00000235371",
    "ENSG00000261083",
    "ENSG00000252999",
    "ENSG00000236576",
    "ENSG00000225518",
    "ENSG00000244402",
    "ENSG00000285571",
    "ENSG00000212374",
    "ENSG00000251324",
    "ENSG00000258824",
    "ENSG00000248508",
    "ENSG00000184574",
    "ENSG00000070831",
    "ENSG00000140254",
    "ENSG00000266985",
    "ENSG00000270741",
    "ENSG00000223390",
    "ENSG00000233862",
    "ENSG00000165233",
    "ENSG00000231889",
    "ENSG00000248428",
    "ENSG00000125868",
    "ENSG00000286180",
    "ENSG00000273084",
    "ENSG00000029364",
    "ENSG00000255474",
    "ENSG00000237256",
]

# The number of genes to use for each query
NUM_GENES = [1, 5, 10, 25, 50, 100, 200]

# The groupby options to use for each query
GROUPBY_OPTIONS = [None, "disease", "publication"]

# Construct the POST bodies for each combination of groupby and number of genes
POST_BODIES = []
for groupby in GROUPBY_OPTIONS:
    for num_genes in NUM_GENES:
        payload = {
            "filter": {
                "dataset_ids": [],
                "development_stage_ontology_term_ids": [],
                "disease_ontology_term_ids": [],
                "gene_ontology_term_ids": RANDOM_GENES[:num_genes],
                "organism_ontology_term_id": "NCBITaxon:9606",
                "self_reported_ethnicity_ontology_term_ids": [],
                "sex_ontology_term_ids": [],
                "publication_citations": [],
            },
            "is_rollup": True,
        }
        if groupby is not None:
            payload["compare"] = groupby

        POST_BODIES.append(payload)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run WMG performance profiling against a deployed environment")
    parser.add_argument("api_url", type=str, help="The WMG API url to query")
    args = parser.parse_args()
    api_url = args.api_url
    profiling_dicts = []
    for i, body in enumerate(POST_BODIES):
        print(f"Executing query {i+1}/{len(POST_BODIES)} with body: {body}")

        response = requests.post(
            f"{api_url}/wmg/v2/query", data=json.dumps(body), headers={"Content-Type": "application/json"}
        )
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
        response_time = response.elapsed.total_seconds()
        download_time = response_time - total_time
        profiling_dict["total_time_backend"] = total_time
        profiling_dict["total_time_response"] = response_time
        profiling_dict["download_time"] = download_time
        response_size = len(response.content) / (1024 * 1024)
        profiling_dict["response_size_mb"] = response_size

        print("Profiling results:", profiling_dict)

        profiling_dict["params"] = {
            "payload": body,
            "num_genes": len(body["filter"]["gene_ontology_term_ids"]),
            "compare": body.get("compare"),
            "api_url": api_url,
        }
        profiling_dicts.append(profiling_dict)

    # Write the profiling results to a file
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    with open(f"profiling_results_{timestamp}.json", "w") as f:
        json.dump(profiling_dicts, f)
