import json


def output_json(data, path):
    with open(path, "w") as f:
        json.dump(data, f)
