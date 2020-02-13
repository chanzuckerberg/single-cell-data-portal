import os, logging

from chalice import Chalice


app = Chalice(app_name="browser")
logger = logging.getLogger(__name__)

swagger_spec_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), f'{os.environ["APP_NAME"]}-api.yml')


@app.route("/")
def index():
    app.log.debug("This is a debug logging test with the app logger")
    app.log.error("This is an error logging test with the app logger")
    logger.debug("This is a debug logging test with the module logger")
    logger.error("This is an error logging test with the module logger")
    return {"hello": "world"}


@app.route("/projects",  methods=['GET'])
def projects():
    return 200


@app.route("/projects/{project_uuid}",  methods=['GET'])
def get_project(project_uuid):
    return 200


@app.route("/projects/{project_uuid}/files/{file_uuid}", methods=['POST'])
def post_file(file_uuid):
    return 200
