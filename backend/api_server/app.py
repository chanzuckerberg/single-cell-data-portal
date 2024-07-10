from backend.common.server.config import create_api_app
from backend.common.utils.json import CurationJSONEncoder
from backend.gene_info.api.ensembl_ids import GeneChecker

app = create_api_app(
    api_paths_and_spec_files=[
        ("/dp", "portal/api/portal-api.yml"),
        ("/cellguide", "cellguide/api/cellguide-api.yml"),
        ("/curation", "curation/api/curation-api.yml"),
        ("/gene_info", "gene_info/api/gene-info-api.yml"),
    ],
    static_folder="backend/api_server/static",
)
app.apis["/curation"].blueprint.json_encoder = CurationJSONEncoder

GeneChecker()

if __name__ == "__main__":
    app.run(host="0.0.0.0", debug=True)
