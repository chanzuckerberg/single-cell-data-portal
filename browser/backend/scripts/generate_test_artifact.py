"""
Generates a stripped-down test JSON artifact to be stored in `tests/unit/dcp-prototype/backend/browser/fixtures`
- Reads an existing JSON artifact from S3
- Selects a subset of files and their respectives projects
- Writes a new `test_artifact.json` to the current directory
"""

import json

import boto3

s3 = boto3.client("s3")
s3.download_file("dcp-test_artifacts", "Artifact.Mar18.json", "artifact.json")

with open("artifact.json", "r") as f:
    content = json.loads(f.read())

    files = content["files"][:3]
    pids = [file["project_id"] for file in files]
    projects = [project for project in content["projects"] if project["id"] in pids]

    content["projects"] = projects
    content["files"] = files

with open("test_artifact.json", "w") as f:
    f.write(json.dumps(content, indent=2))
