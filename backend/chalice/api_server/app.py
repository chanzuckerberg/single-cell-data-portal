import os
import sys
from functools import wraps

import chalice
from chalice import Chalice, CORSConfig

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.utils.db_utils import DbUtils
from backend.corpora.common.utils.s3_utils import generate_file_url
from backend.corpora.common.authorizer import assert_authorized

app = Chalice(app_name=f"{os.environ['APP_NAME']}-{os.environ['DEPLOYMENT_STAGE']}")

cors_config = CORSConfig(allow_origin="*", max_age=600, allow_credentials=True)


def requires_auth():
    """
    A decorator for assert_authorized
    :return: 401 or the original function response.
    """

    def decorate(func):
        @wraps(func)
        def call(*args, **kwargs):
            assert_authorized(app.current_request.headers)
            return func(*args, **kwargs)

        return call

    return decorate


@app.route("/")
def index():
    return "render documentation here"


@app.route("/projects", cors=cors_config)
def get_projects():
    db = DbUtils()
    projects = db.query_projects()

    return chalice.Response(status_code=200, headers={"Content-Type": "application/json"}, body=projects)


@app.route("/projects/{project_id}", cors=cors_config)
def get_project(project_id: str):
    db = DbUtils()
    project = db.query_project(project_id)

    return chalice.Response(
        status_code=200 if project else 404, headers={"Content-Type": "application/json"}, body=project,
    )


@app.route("/projects/{project_id}/files", cors=cors_config)
@requires_auth()
def get_project_files(project_id: str):
    db = DbUtils()
    files = db.query_downloadable_project_files(project_id)
    return chalice.Response(
        status_code=200 if files else 404, headers={"Content-Type": "application/json"}, body=files,
    )


@app.route("/files/{file_id}", cors=cors_config)
@requires_auth()
def get_file(file_id: str):
    db = DbUtils()
    file = db.query_file(file_id)

    download_url = ""
    if file:
        project = db.query_project(file["project_id"])
        file_prefix = f"{project['label']}/matrix.loom"
        download_url = generate_file_url(file_prefix)

    return chalice.Response(
        status_code=200 if file else 404, headers={"Content-Type": "application/json"}, body={"url": download_url},
    )
