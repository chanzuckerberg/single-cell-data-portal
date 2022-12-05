from os import popen

PROJECT_ROOT = popen("git rev-parse --show-toplevel").read().strip()
if PROJECT_ROOT == "":
    FIXTURES_ROOT = "tests/unit/backend/wmg/fixtures"
else:
    FIXTURES_ROOT = PROJECT_ROOT + "/tests/unit/backend/wmg/fixtures"
