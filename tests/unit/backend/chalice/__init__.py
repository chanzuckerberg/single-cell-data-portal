import functools

import requests
from chalice.cli import CLIFactory
from chalice.local import LocalGateway, LocalGatewayException

import subprocess


def run(args):
    print(f"RUN: {' '.join(args)}")
    response = subprocess.run(args, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if response.stdout:
        print(f'RUN STDOUT:\n{response.stdout.decode("utf-8")}')
    if response.stderr:
        print(f'RUN STDERR:\n{response.stderr.decode("utf-8")}')
    print("\n")
    return response.stdout.decode("utf-8")


class ChaliceTestHarness:
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
