"""
After a PR is merged, move the corresponding issue from the "Code Review" pipeline to the "Ready for Prod" pipeline.
"""

import os
from typing import List

import requests

access_token = os.environ["zenhub_api_key"]
endpoint = "https://api.zenhub.com/public/graphql"
headers = {"Authorization": f"Bearer {access_token}"}
pipeline_id = "Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzI5NzY2ODc"  # "Code Review" pipeline in "single-cell" workspace
repo_id = "Z2lkOi8vcmFwdG9yL1JlcG9zaXRvcnkvNTM5NzQwOTg"  # "single-cell-data-portal" repo
destination_pipeline = "Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzIyMDg1MDE"  # "Ready for Prod" pipeline in "single-cell" workspace


def move_to_ready_for_prod() -> None:
    # Get all issues in the Code Review pipeline.
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
    [issue["id"] for issue in issues]
    moving_issues: List[str] = []
    moving_issue_strings: List[str] = []
    for issue in issues:
        if issue["connectedPrs"]["nodes"]:
            # Filter out issue if it has any open PRs connected to it.
            open_prs = [
                f"PR #{i['number']}: {i['title']}" for i in issue["connectedPrs"]["nodes"] if i["state"] == "OPEN"
            ]
            if open_prs:
                print(f"Issue #{issue['number']} has open PRs:", *open_prs, sep="\n\t")
                continue
        moving_issues.append(issue["id"])
        moving_issue_strings.append(f"issue #{issue['number']}: {issue['title']}")

    # Move issues.
    print("Moving issues:", *moving_issue_strings, sep="\n\t")
    mutate_move_isses = """mutation movePipelineIssues {
      movePipelineIssues(input: {
          pipelineId: "Z2lkOi8vcmFwdG9yL1BpcGVsaW5lLzIyMDg1MDE"
          issueIds: $MOVING_ISSUES
      }) {
          successCount
      }
    }""".replace(
        "$MOVING_ISSUES", str(moving_issues)
    )
    resp = requests.post(endpoint, json={"query": mutate_move_isses}, headers=headers)
    resp.raise_for_status()
    print(resp.json())


if __name__ == "__main__":
    move_to_ready_for_prod()
