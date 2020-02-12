#!/usr/bin/env python

import os, sys, json, unittest, functools
import requests
from chalice.cli import CLIFactory
from chalice.local import LocalGateway, LocalGatewayException
from furl import furl


class ChaliceTestHarness:
    def __init__(self):
        project_dir = os.path.join(os.path.dirname(__file__), "..")
        config = CLIFactory(project_dir=project_dir).create_config_obj(chalice_stage_name="dev")
        self._chalice_app = config.chalice_app
        self._gateway = LocalGateway(self._chalice_app, config)

    @functools.lru_cache(maxsize=64, typed=False)
    def __getattr__(self, method):
        return functools.partial(self.request, method=method.upper())

    def request(self, path, headers={}, body={}, method="GET"):
        resp_obj = requests.Response()
        try:
            response = self._gateway.handle_request(method, path, headers, body)
        except LocalGatewayException as error:
            resp_obj.status_code = error.CODE
            resp_obj.headers = error.headers
            resp_obj.body = error.body
        else:
            resp_obj.status_code = response['statusCode']
            resp_obj.headers = response['headers']
            resp_obj.body = response['body']
        resp_obj.headers['Content-Length'] = str(len(body))
        return resp_obj


class TestChaliceApp(unittest.TestCase):
    def setUp(self):
        self.app = ChaliceTestHarness()

    def test_root_route(self):
        res = self.app.get('/')
        res.raise_for_status()
        self.assertEqual(res.status_code, requests.codes.ok)


if __name__ == '__main__':
    unittest.main()
