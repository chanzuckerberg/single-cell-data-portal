"""
This script closes issues in the Ready for Prod pipeline that are not blocked by open issues or PRs.
"""
import os
from typing import List

import requests

access_token = os.environ["zenhub_api_key"]
endpoint = "https://api.zenhub.com/public/graphql"
headers = {"Authorization": f"Bearer {access_token}"}
pipeline_id = "Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzIyMDg1MDE"  # "Ready for Prod" pipeline in "single-cell" workspace
repo_id = "Z2lkOi8vcmFwdG9yL1JlcG9zaXRvcnkvNTM5NzQwOTg"  # "single-cell-data-portal" repo


def close_ready_for_prod() -> None:
    # Get all issues in the Ready for Prod pipeline.
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
    issues = resp.json()["data"]["searchIssuesByPipeline"]["nodes"]

    # Filter out issues that are blocked by open issues or PRs.
    issue_ids = [issue["id"] for issue in issues]
    closing_issues: List[str] = []
    closing_issue_strings: List[str] = []
    for issue in issues:
        if issue["connectedPrs"]["nodes"]:
            # Filter out issue if it has any open PRs connected to it.
            open_prs = [
                f"PR #{i['number']}: {i['title']}" for i in issue["connectedPrs"]["nodes"] if i["state"] == "OPEN"
            ]
            if open_prs:
                print(f"issue #{issue['number']} has open PRs:", *open_prs, sep="\n\t")
                continue
        if issue["blockingIssues"]["nodes"]:
            # Filter out issue if it is blocked by any open issues this is not also in the ready for prod pipeline.
            blocking_issues = [
                f"issue #{i['number']}: {i['title']}"
                for i in issue["blockingIssues"]["nodes"]
                if i["state"] == "OPEN" and i["id"] != issue_ids
            ]
            if blocking_issues:
                print(f"Issue {issue['number']} is blocked by:", *blocking_issues, sep="\n\t")
                continue
            closing_issues.append(issue["id"])
            closing_issue_strings.append(f"issue #{issue['number']}: {issue['title']}")

    # Close issues.
    print("Closing issues:", *closing_issue_strings, sep="\n\t")
    mutate_close_isses = """mutation closeIssues {
      closeIssues(input: {
          issueIds: $CLOSING_ISSUES
      }) {
          successCount
      }
    }""".replace(
        "$CLOSING_ISSUES", str(closing_issues)
    )
    resp = requests.post(endpoint, json={"query": mutate_close_isses}, headers=headers)
    resp.raise_for_status()
    print(resp.json())


if __name__ == "__main__":
    close_ready_for_prod()
