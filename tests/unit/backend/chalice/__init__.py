import functools

import requests
from chalice.cli import CLIFactory
from chalice.local import LocalGateway, LocalGatewayException


class ChaliceTestHarness:
    """
    This class initializes a chalice application that can be used to make API requests against the application.
    """

    def __init__(self, project_dir):
        config = CLIFactory(project_dir=project_dir).create_config_obj(chalice_stage_name="dev")
        self._chalice_app = config.chalice_app
        self._gateway = LocalGateway(self._chalice_app, config)

    @functools.lru_cache(maxsize=128, typed=False)
    def __getattr__(self, item):
        item = item.upper()
        return functools.partial(self.request, method=item)

    def request(self, path, headers=None, data="", method="GET"):
        resp_obj = requests.Response()
        if not headers:
            headers = {}
        try:
            response = self._gateway.handle_request(method, path, headers, data)
        except LocalGatewayException as error:
            resp_obj.status_code = error.CODE
            resp_obj.headers = error.headers
            resp_obj.body = error.body
        else:
            resp_obj.status_code = response["statusCode"]
            resp_obj.headers = response["headers"]
            resp_obj.body = response["body"]
        resp_obj.headers["Content-Length"] = str(len(data))
        return resp_obj
