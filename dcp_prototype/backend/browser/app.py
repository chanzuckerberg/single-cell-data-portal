import logging

from chalice import Chalice

app = Chalice(app_name="browser")
logger = logging.getLogger(__name__)


@app.route("/")
def index():
    app.log.debug("This is a debug logging test with the app logger")
    app.log.error("This is an error logging test with the app logger")
    logger.debug("This is a debug logging test with the module logger")
    logger.error("This is an error logging test with the module logger")
    return {"hello": "world"}


@app.route("/projects", methods=["GET"])
def projects():
    return {"projects": []}


@app.route("/projects/{project_uuid}", methods=["GET"])
def get_project(project_uuid):
    return {"project_uuid": {}}


@app.route("/files/{file_uuid}", methods=["POST"])
def post_file(file_uuid):
    return {"hello": "world"}
