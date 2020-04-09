import os
import sys

import chalice
from chalice import Chalice

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.code.common.db_utils import DbUtils
from browser.code.common.s3_utils import generate_file_url

app = Chalice(app_name=f"{os.environ['APP_NAME']}-{os.environ['DEPLOYMENT_STAGE']}")


@app.route("/")
def index():
    return "render documentation here"


@app.route("/projects")
def get_projects():
    db = DbUtils()
    projects = db.query_projects()

    return chalice.Response(
        status_code=200, headers={"Content-Type": "application/json", "Access-Control-Allow-Origin": "*"}, body=projects
    )


@app.route("/projects/{project_id}")
def get_project(project_id):
    db = DbUtils()
    project = db.query_project(project_id)

    return chalice.Response(
        status_code=200 if project else 404,
        headers={"Content-Type": "application/json", "Access-Control-Allow-Origin": "*"},
        body=project,
    )


@app.route("/projects/{project_id}/files")
def get_project_files(project_id):
    db = DbUtils()
    files = db.query_downloadable_project_files(project_id)

    return chalice.Response(
        status_code=200 if files else 404,
        headers={"Content-Type": "application/json", "Access-Control-Allow-Origin": "*"},
        body=files,
    )


@app.route("/files/{file_id}")
def get_file(file_id):
    db = DbUtils()
    file = db.query_file(file_id)

    download_url = ""
    if file:
        project = db.query_project(file.project_id)
        file_prefix = f"{project['label']}/matrix.loom"
        download_url = generate_file_url(file_prefix)

    return chalice.Response(
        status_code=200 if file else 404,
        headers={"Content-Type": "application/json", "Access-Control-Allow-Origin": "*"},
        body={"url": download_url},
    )
