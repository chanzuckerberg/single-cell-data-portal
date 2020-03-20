import os
import sys

import chalice
from chalice import Chalice

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.code.common.browser_orm import DBSessionMaker, Project, File
from browser.code.common.db_utils import (
    get_project_assays,
    get_project_organs,
    get_project_species,
    get_project_contributors,
    get_downloadable_project_files,
)
from browser.code.common.s3_utils import generate_file_url

app = Chalice(app_name=f"{os.environ['APP_NAME']}-{os.environ['DEPLOYMENT_STAGE']}")
session = DBSessionMaker().session()


@app.route("/")
def index():
    return "render documentation here"


@app.route("/projects")
def get_projects():
    projects = []
    for project in session.query(Project).all():
        projects.append(
            {
                "id": project.id,
                "title": project.title,
                "assays": get_project_assays(project.id),
                "organs": get_project_organs(project.id),
                "species": get_project_species(project.id),
                "cell_count": project.cell_count,
            }
        )

    return chalice.Response(
        status_code=200, headers={"Content-Type": "application/json", "Access-Control-Allow-Origin": "*"}, body=projects
    )


@app.route("/projects/{project_id}")
def get_project(project_id):
    project = session.query(Project).get(project_id)

    response_body = {
        "id": project.id,
        "title": project.title,
        "assays": get_project_assays(project.id),
        "organs": get_project_organs(project.id),
        "species": get_project_species(project.id),
        "contributors": get_project_contributors(project.id),
        "description": project.description,
        "biosample_categories": project.biosample_categories.split(","),
        "development_stages": project.development_stages.split(","),
        "diseases": project.diseases.split(","),
        "cell_isolation_methods": project.cell_isolation_methods.split(","),
        "cell_types": project.cell_types.split(","),
        "cell_count": project.cell_count,
        "paired_end": project.paired_end.split(","),
        "nucleic_acid_sources": project.nucleic_acid_sources.split(","),
        "input_nucleic_acid_molecules": project.input_nucleic_acid_molecules.split(","),
        "publication_title": project.publication_title,
        "publication_doi": project.publication_doi,
    }

    return chalice.Response(
        status_code=200,
        headers={"Content-Type": "application/json", "Access-Control-Allow-Origin": "*"},
        body=response_body,
    )


@app.route("/projects/{project_id}/files")
def get_project_files(project_id):
    files = get_downloadable_project_files(project_id, session)

    return chalice.Response(
        status_code=200, headers={"Content-Type": "application/json", "Access-Control-Allow-Origin": "*"}, body=files
    )


@app.route("/files/{file_id}")
def get_file(file_id):
    file = session.query(File).get(file_id)
    project = session.query(Project).get(file.project_id)

    # file_prefix = f"{project.title}/{file.filename}"
    file_prefix = f"{project.title}/matrix.loom"
    download_url = generate_file_url(file_prefix)

    response_body = {"url": download_url}
    return chalice.Response(
        status_code=200,
        headers={"Content-Type": "application/json", "Access-Control-Allow-Origin": "*"},
        body=response_body,
    )
