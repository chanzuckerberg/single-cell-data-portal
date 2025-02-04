import json
import os

import jsonschema

local_path = os.path.dirname(os.path.realpath(__file__))
with open(local_path + "/schema.json") as fp:
    schema = json.load(fp)
validator = jsonschema.Draft7Validator(schema)
