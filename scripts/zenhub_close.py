"""
This script closes issues in the Ready for Prod pipeline that are not blocked by open issues or PRs.
"""
import json
import logging
import os
import sys
from typing import Iterable, List, Tuple

import requests

logging.basicConfig(level=logging.INFO, format="%(message)s", handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger(__name__)

pipeline_id = "Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzMwNTE5ODA"  # "Done" pipeline in "single-cell" workspace
repo_id = "Z2lkOi8vcmFwdG9yL1JlcG9zaXRvcnkvMTMzNDEyNzQ1"  # "Zenhub-experiment" repo


def get_issues(pipeline_id: str, repo_id: str) -> List[dict]:
    query_list_issues = """query {
    searchIssuesByPipeline(
        pipelineId: "$PIPELINE_ID",
        filters: {
          repositoryIds: "$REPO_ID"
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
        "$REPO_ID", repo_id
    )
    resp = requests.post(endpoint, json={"query": query_list_issues}, headers=headers)
    resp.raise_for_status()
    return resp.json()["data"]["searchIssuesByPipeline"]["nodes"]


def close_issues(issue_ids: Iterable[str]) -> dict:
    mutate_close_isses = """mutation closeIssues {
      closeIssues(input: {
          issueIds: $CLOSING_ISSUES
      }) {
          successCount
      }
    }""".replace(
        "$CLOSING_ISSUES", json.dumps(issue_ids)
    )
    resp = requests.post(endpoint, json={"query": mutate_close_isses}, headers=headers)
    resp.raise_for_status()
    return resp.json()


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


def close_ready_for_prod() -> None:
    # Get all issues in the Ready for Prod pipeline.
    issues = get_issues(pipeline_id, repo_id)

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
        logger.info("Closing issues:\n\t" + "\n\t".join(issue_strings))
        data = close_issues(issue_ids)
        logger.info(data)
        assert data["data"]["closeIssues"]["successCount"] == len(issues_to_close)


if __name__ == "__main__":
    access_token = os.environ["ZENHUB_TOKEN"]
    endpoint = "https://api.zenhub.com/public/graphql"
    headers = {"Authorization": f"Bearer {access_token}"}
    close_ready_for_prod()
