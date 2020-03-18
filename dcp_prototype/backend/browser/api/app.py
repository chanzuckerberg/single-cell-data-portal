import os
import sys

import chalice
from chalice import Chalice

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.code.common.browser_orm import DBSessionMaker, Project, File
from browser.code.common.db_utils import (
    get_project_assays,
    get_project_tissues,
    get_project_species,
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
                "tissues": get_project_tissues(project.id),
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
        "label": project.label,
        "assays": get_project_assays(project.id),
        "tissues": get_project_tissues(project.id),
        "species": get_project_species(project.id),
        "description": project.description,
        "category": project.category,
        "developmental_stage": project.developmental_stage,
        "disease_ontology": project.disease_ontology,
        "sample_type": project.sample_type,
        "organ_part": project.organ_part,
        "analysis_protocol": project.analysis_protocol.split(","),
        "cell_count": project.cell_count,
        "donor_count": project.donor_count,
        "publication_title": project.publication_title,
        "publication_doi": project.publication_doi,
        "contact_name": project.contact_name,
        "contact_institution": project.contact_institution,
        "contact_email": project.contact_email,
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

    file_prefix = f"{file.project_id}/{file.filename}"
    download_url = generate_file_url(file_prefix)

    response_body = {"url": download_url}
    return chalice.Response(
        status_code=200,
        headers={"Content-Type": "application/json", "Access-Control-Allow-Origin": "*"},
        body=response_body,
    )
