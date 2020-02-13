import os, logging

from chalice import Chalice

log_level = logging.getLevelName(os.environ.get("LOG_LEVEL", "INFO"))
logging.basicConfig(level=log_level)
logging.getLogger().setLevel(log_level)

app = Chalice(app_name="browser")
app.debug = True if log_level == logging.DEBUG else False
for l in (
    "botocore.vendored.requests.packages.urllib3.connectionpool",
    "requests.packages.urllib3.connectionpool",
):
    logging.getLogger(l).setLevel(max(log_level, logging.WARNING))

logger = logging.getLogger(__name__)


@app.route("/")
def index():
    app.log.debug("This is a debug logging test with the app logger")
    app.log.error("This is an error logging test with the app logger")
    logger.debug("This is a debug logging test with the module logger")
    logger.error("This is an error logging test with the module logger")
    return {"hello": "world"}
