#!/usr/bin/env python

import json
import requests


class WarmCache():
    def __init__(self):
        self.dev_url_base = "https://api.cellxgene.dev.single-cell.czi.technology/"
        self.staging_url_base = "https://api.cellxgene.staging.single-cell.czi.technology/"
        self.prod_url_base = "https://api.cellxgene.cziscience.com/"

        self.dev_project_ids = self.get_project_ids(self.dev_url_base)
        self.staging_project_ids = self.get_project_ids(self.staging_url_base)
        self.prod_project_ids = self.get_project_ids(self.prod_url_base)

    def get_project_ids(self, url_base):
        response = requests.get(f"{url_base}dp/v1/project")
        project_list = json.loads(response._content)
        return [x["id"] for x in project_list["projects"]]

    def get_project_info(self, url_base, project_id):
        response = requests.get(f"{url_base}dp/v1/project/{project_id}")

    def get_all_projects(self):
        for project_id in self.dev_project_ids:
            self.get_project_info(self.dev_url_base, project_id)

        for project_id in self.staging_project_ids:
            self.get_project_info(self.staging_url_base, project_id)

        for project_id in self.prod_project_ids:
            self.get_project_info(self.prod_url_base, project_id)


if __name__ == "__main__":
    WarmCache().get_all_projects()
