"""
A simple script demonstrating how to access the deployed database using the corpora application layer.

Follow the Development instructions up to step 4 before running this script.
https://github.com/chanzuckerberg/corpora-data-portal/tree/main/backend/chalice/api_server#development
"""

import os
import sys
from pprint import pprint

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.entities import Project, Dataset

os.environ["DEPLOYMENT_STAGE"] = 'dev'
os.environ["CORPORA_LOCAL_DEV"] = ''


projects = Project.list()
print(projects)
pprint(projects[0].to_dict())

