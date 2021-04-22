#!/usr/bin/env python3
import dateutil.parser as dp
import json
import os
import sys

import click
import requests

github_org = "chanzuckerberg"
github_repo = "corpora-data-portal"
github_graphql_endpoint = "https://api.github.com/graphql"
github_deployment_endpoint = "https://api.github.com/repos/chanzuckerberg/corpora-data-portal/deployments"


def get_latest_successful_deployment(github_api_token, stage):
    """get the latest successful/active deployment github sha"""
    # Assumption: One of the most recent 50 deployment attempts was successful
    query = """
    query($repo_owner:String!, $repo_name:String!, $deployment_env:String!) {
        repository(owner: $repo_owner, name: $repo_name) {
          deployments(environments: [$deployment_env], last: 50) {
            nodes {
              commitOid
              statuses(first: 100) {
                nodes {
                  state
                  updatedAt
                }
              }
            }
          }
        }
      }
      """

    variables = {"repo_owner": github_org, "repo_name": github_repo, "deployment_env": stage}

    headers = {"Authorization": "token %s" % github_api_token}
    query = {"query": query, "variables": variables}

    try:
        resp = requests.post(url=github_graphql_endpoint, json=query, headers=headers)
        if resp.status_code != 200:
            print("Error: Unexpected response {}".format(resp))
            return None
    except requests.exceptions.RequestException as e:
        print("Error: {}".format(e))
        return None

    resp_json = json.loads(resp.text)
    deployments = resp_json["data"]["repository"]["deployments"]

    sha_tuple = (None, None)

    for node in deployments["nodes"]:
        gh_sha = node["commitOid"]
        for status in node["statuses"]["nodes"]:
            if status["state"] == "SUCCESS":
                parsed_t = dp.parse(status["updatedAt"])
                if sha_tuple[0] == None:
                    sha_tuple = (gh_sha, parsed_t)
                else:
                    if sha_tuple[1] < parsed_t:
                        sha_tuple = (gh_sha, parsed_t)
                break

    return sha_tuple


def trigger_deploy(github_api_token, deployment_stage, github_sha, dry_run):
    """Start deployment to the given environment based on the github sha"""
    headers = {"Authorization": "token %s" % github_api_token, "Accept": "application/vnd.github.v3.text-match+json"}

    tag = f"sha-{github_sha[0:8]}"

    params = {
        "ref": github_sha,
        "auto_merge": False,
        "environment": deployment_stage,
        "required_contexts": [],
        "payload": {"tag": tag},
    }

    if dry_run:
       print(f"Dry run requested. Would deploy {tag} to environment {deployment_stage}")
       return

    print(f"Deploying {tag} to environment {deployment_stage}")
    try:
        resp = requests.post(github_deployment_endpoint, headers=headers, json=params)
        if resp.status_code != 201:
            print("Error: Unexpected response {}".format(resp))
            return
    except requests.exceptions.RequestException as e:
        print("Error: {}".format(e))
        return

    print("Deployment successful")


@click.command()
@click.argument("deployment_stage")
@click.option("--github_sha", help="github sha to be deployed", default=None)
@click.option("--dry_run", help="do not perform actual deployment", default=False, is_flag=True)
def happy_deploy(deployment_stage, github_sha, dry_run):
    api_token = os.getenv("GITHUB_TOKEN")
    if api_token is None:
        print("Error: Please set GITHUB_TOKEN environment variable")
        return

    deployment_stage = "stage"

    # If github sha is not provided, get the latest succesful deployment
    # github sha of staging environment
    if github_sha is None:
        github_sha, parsed_t = get_latest_successful_deployment(api_token, deployment_stage)
        print(f"Latest succesful '{deployment_stage}' deployment on {parsed_t}: commit {github_sha}")

    if github_sha is None:
        print(f"Error: Could not find a successful deployment for deployment stage {deployment_stage}, and no --github_sha was given")
        sys.exit(1)

    # Trigger deployment on the given stage. This will trigger github actions
    # and start/update the deployment.
    if github_sha is not None:
        trigger_deploy(api_token, deployment_stage, github_sha, dry_run)


if __name__ == "__main__":
    happy_deploy()
