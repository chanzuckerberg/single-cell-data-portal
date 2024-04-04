"""
This script closes issues in the Ready for Prod pipeline that are not blocked by open issues or PRs.

Use zenhubs explorer to test https://developers.zenhub.com/explorer
"""
import json
import logging
import os
import sys
from typing import Any, Dict, Iterable, List, Optional, Tuple

import requests

logger = logging.getLogger(__name__)

ACCESS_TOKEN = "zh_6dfcfccc7b6668a5463b74661e31f3dcc100a30f97e51801528bf2ca162ab26f"
WORKSPACE_NAME = "single-cell"  # workspace name
REPO_NAMES = ["single-cell-data-portal"]  # list the repos you want to close issues in
SOURCE_PIPELINE_NAME = "Ready for Prod"  # the pipelines you want to close issues in


class ZenHubProvider:
    def __init__(self, endpoint: str = "https://api.zenhub.com/graphql", access_token: Optional[str] = None):
        self.access_token = access_token or os.environ.get("ZENHUB_TOKEN", "")
        self.headers = {"Authorization": f"Bearer {access_token}"}
        self.endpoint = endpoint
        self.session = requests.Session()
        self.session.headers.update(self.headers)

    def get_workspaces(self, workspace_name: str) -> List[Any]:
        """
        Get all workspaces and their repositories and pipelines.
        :param workspace_name: This is required to get all of the workspaces. It should match an existing workspace name.
        :return: a list of workspaces with their repositories and pipelines.
        """
        query_get_id = """query {
      viewer {
        id
        searchWorkspaces(query: "$WORKSPACE_NAME") {
          nodes {
            id
            name
            repositoriesConnection {
              nodes {
                id
                name
              }
            }
            pipelinesConnection{
              nodes{
                id
                name
              }
            }
          }
        }
      }
    }""".replace(
            "$WORKSPACE_NAME", workspace_name
        )
        return self._query(query_get_id)["data"]["viewer"]["searchWorkspaces"]["nodes"]

    def get_issues(self, pipeline_id: str, repo_ids: List[str]) -> List[dict]:
        """
        Get all issues in the Ready for Prod pipeline from the repositories listed
        :param pipeline_id: the pipeline id to search for issues in
        :param repo_ids: the repository ids to search for issues in
        :return: issues in the Ready for Prod pipeline from the repositories listed
        """
        query_list_issues = """query {
        searchIssuesByPipeline(
            pipelineId: "$PIPELINE_ID",
            filters: {
              repositoryIds: $REPO_IDS
            }
        ) {
          nodes {
              id
              number
              title
              blockingIssues{
                nodes{
                  id
                  number
                  title 
                  state
                }
              }
              connectedPrs{
                nodes{
                  id
                  number
                  title
                  state
                }
              }
          }
        }
        }""".replace(
            "$PIPELINE_ID", pipeline_id
        ).replace(
            "$REPO_ID", json.dumps(repo_ids)
        )
        return self._query(query_list_issues)["data"]["searchIssuesByPipeline"]["nodes"]

    def close_issues(self, issue_ids: Iterable[str]) -> dict:
        """
        Close issues by their ids.
        :param issue_ids: the ids of the issues to close
        :return: the response from the server
        """
        mutate_close_isses = """mutation closeIssues {
          closeIssues(input: {
              issueIds: $CLOSING_ISSUES
          }) {
              successCount
          }
        }""".replace(
            "$CLOSING_ISSUES", json.dumps(issue_ids)
        )
        return self._query(mutate_close_isses)

    def _query(self, query: str) -> dict:
        resp = self.session.post(self.endpoint, json={"query": query})
        resp.raise_for_status()
        return resp.json()


def parse_workspace(workspaces: List[Any], workspace_name: str) -> Dict[str, Any]:
    """
    Parse the workspace from the list of nodes.
    :param workspaces: nodes from the response
    :param workspace_name: the workspace name to search for
    :return: the workspace
    """
    for node in workspaces:
        if node["name"] == workspace_name:
            return node
    else:
        raise ValueError(f"Workspace {workspace_name} not found.")


def parse_repo_ids(workspace: Dict[str, Any], repo_names: List[str]) -> List[str]:
    """
    Parse the repository ids from the workspace.
    :param workspace: the workspace_name
    :param repo_names:
    :return:
    """
    repo_ids = []
    for repo in workspace["repositoriesConnection"]["nodes"]:
        if repo["name"] in repo_names:
            repo_ids.append(repo["id"])
    if len(repo_ids) != len(repo_names):
        raise ValueError(f"Only found {len(repo_ids)} out of {len(repo_names)} repos.")
    return repo_ids


def parse_pipeline(workspace: Dict[str, Any], source_pipeline: str) -> str:
    for pipeline in workspace["pipelinesConnection"]["nodes"]:
        if pipeline["name"] == source_pipeline:
            return pipeline["id"]
    else:
        raise ValueError(f"Pipeline {source_pipeline} not found.")


def format_issue(issue: dict) -> str:
    return f"#{issue['number']} ({issue['title']})"


def filter_issues(issues: List[dict]) -> Tuple[List[Tuple[str, str]], List[dict]]:
    # Filter out issues that are blocked by open issues or PRs.
    issue_ids = [issue["id"] for issue in issues]
    issues_to_close: List[Tuple[str, str]] = []
    blocked_issues: List[dict] = []
    for issue in issues:
        open_prs = []
        blocking_issues = []
        issue_string = format_issue(issue)
        if issue["connectedPrs"]["nodes"]:
            # Filter out issue if it has any open PRs connected to it.
            open_prs = [format_issue(i) for i in issue["connectedPrs"]["nodes"] if i["state"] == "OPEN"]
        if issue["blockingIssues"]["nodes"]:
            # Filter out issue if it is blocked by any open issues this is not also in the ready for prod pipeline.
            blocking_issues = [
                format_issue(i)
                for i in issue["blockingIssues"]["nodes"]
                if i["state"] == "OPEN" and i["id"] not in issue_ids
            ]
        if blocking_issues or open_prs:
            blocked_issues.append({"issue": issue_string, "open_prs": open_prs, "blocking_issues": blocking_issues})
        else:
            issues_to_close.append((issue["id"], issue_string))
    return issues_to_close, blocked_issues


def close_ready_for_prod(
    workspace_name: str, repo_names: List[str], source_pipeline_name: str, provider: Optional[ZenHubProvider] = None
) -> None:
    provider = provider or ZenHubProvider()

    # translate names into zenhub ids
    workspace_resp = provider.get_workspaces(workspace_name)
    workspace = parse_workspace(workspace_resp, workspace_name)
    repo_ids = parse_repo_ids(workspace, repo_names)
    pipeline_id = parse_pipeline(workspace, source_pipeline_name)

    # Get all issues in the Ready for Prod pipeline.
    issues = provider.get_issues(pipeline_id, repo_ids)

    # Filter out issues that are blocked by open issues or PRs.
    issues_to_close, blocked_issues = filter_issues(issues)

    # Blocked issues.
    if blocked_issues:
        for blocked_issue in blocked_issues:
            logger.info(
                f"Not closed: {blocked_issue['issue']}"
                + "\n\tBlocked by:\n\t\t"
                + "\n\t\t".join(blocked_issue["blocking_issues"])
                + "\n\tOpen PRs:\n\t\t"
                + "\n\t\t".join(blocked_issue["open_prs"])
            )

    # Close issues.
    if not issues_to_close:
        logger.info("No issues to close.")
    else:
        issue_ids, issue_strings = zip(*issues_to_close)
        logger.info("Closing issues:\n\t" + "- ".join(issue_strings))
        data = provider.close_issues(issue_ids)
        logger.info(data)
        assert data["data"]["closeIssues"]["successCount"] == len(issues_to_close)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=[logging.StreamHandler(sys.stdout)])
    try:
        close_ready_for_prod(WORKSPACE_NAME, REPO_NAMES, SOURCE_PIPELINE_NAME)
    except Exception as e:
        logger.exception(e)
        sys.exit(1)
