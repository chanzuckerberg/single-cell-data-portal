import concurrent.futures
import sys

import dcplib.etl
import hca


DSS_SEARCH_QUERY_TEMPLATE = {
    "query": {
        "bool": {
            "must": [
                {
                    "bool": {
                        "should": [],
                        "minimum_number_should_match": 1
                    }
                },
                {
                    "bool": {
                        "should": [
                            {
                                "match": {
                                    "files.donor_organism_json.biomaterial_core.ncbi_taxon_id": 9606
                                }
                            },
                            {
                                "match": {
                                    "files.donor_organism_json.biomaterial_core.ncbi_taxon_id": 10090
                                }
                            }
                        ],
                        "minimum_number_should_match": 1
                    }
                },
                {
                    "bool": {
                        "should": [
                            {
                                "match": {
                                    "files.analysis_process_json.type.text": "analysis"
                                }
                            },
                            {
                                "match": {
                                    "files.analysis_process_json.process_type.text": "analysis"
                                }
                            }
                        ],
                        "minimum_number_should_match": 1
                    }
                },
            ]
        }
    }
}

project_uuid = sys.argv[1]

q = DSS_SEARCH_QUERY_TEMPLATE
q['query']['bool']['must'][0]['bool']['should'].append({
    "match": {
	"files.project_json.provenance.document_id": project_uuid
    }
})

CLIENT = hca.dss.DSSClient()
extractor = dcplib.etl.DSSExtractor(
    staging_directory='project_data',
    content_type_patterns=['application/json; dcp-type="metadata*"'],
    filename_patterns=["*.fastq*", "*.bam", "*zarr*", "genes.tsv", "barcodes.tsv"],
    dss_client=CLIENT)

extractor.extract(
    query=q,
    max_workers=2,
    dispatch_executor_class=concurrent.futures.ThreadPoolExecutor)
