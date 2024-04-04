from typing import Tuple
from unittest.mock import Mock, patch

import pytest
import requests
from requests import Response

from scripts.zenhub_close import (
    ZenHubProvider,
    close_ready_for_prod,
    filter_issues,
    format_issue,
    parse_pipeline,
    parse_repo_ids,
    parse_workspace,
)


@pytest.fixture
def sample_issues() -> list:
    return [
        {
            "id": 1,
            "number": 101,
            "title": "Sample Issue 1",
            "connectedPrs": {"nodes": []},
            "blockingIssues": {"nodes": []},
            "state": "OPEN",
        }
    ]


@pytest.fixture
def blocking_issue() -> dict:
    return {"id": 3, "state": "OPEN", "number": 301, "title": "Blocked Issue"}


@pytest.fixture
def mock_zenhub_provider_query() -> Tuple[ZenHubProvider, Mock]:
    # Create a mock ZenHubProvider instance with a mocked _query method
    with patch.object(ZenHubProvider, "_query") as mock_query:
        zenhub_provider = ZenHubProvider()
        yield zenhub_provider, mock_query


@pytest.fixture
def mock_workspace_node() -> dict:
    return {
        "id": "workspace_id",
        "name": "example_workspace",
        "repositoriesConnection": {"nodes": [{"id": "repo_id", "name": "repo_name"}]},
        "pipelinesConnection": {"nodes": [{"id": "pipeline_id", "name": "pipeline_name"}]},
    }


@pytest.fixture
def mock_zenhub_provider(mock_workspace_node) -> Mock:
    # Create a mock ZenHubProvider instance with a mocked _query method
    mock = Mock()
    mock.get_workspaces.return_value = [mock_workspace_node]
    mock.get_issues = Mock(return_value=[])
    mock.close_issues = Mock(return_value={"data": {"closeIssues": {"successCount": 0}}})
    return mock


class TestZenHubProvider:
    def test_get_workspace(self, mock_zenhub_provider_query, mock_workspace_node):
        # Prepare mock data and expected response
        workspace_name = "example_workspace"
        _, mock_query = mock_zenhub_provider_query
        mock_query.return_value = {"data": {"viewer": {"searchWorkspaces": {"nodes": [mock_workspace_node]}}}}

        # Call the method
        zenhub_provider, mock_query = mock_zenhub_provider_query
        result = zenhub_provider.get_workspaces(workspace_name)

        # Assert the expected response
        assert result == [mock_workspace_node]
        assert '"example_workspace"' in mock_query.call_args_list[0][0][0]

    def test_get_issues(self, mock_zenhub_provider_query):
        # Prepare mock data and expected response
        pipeline_id = "pipeline_id"
        repo_ids = ["repo_id"]
        mock_response = {
            "data": {"searchIssuesByPipeline": {"nodes": [{"id": "issue_id", "number": 1, "title": "Issue Title"}]}}
        }
        _, mock_query = mock_zenhub_provider_query
        mock_query.return_value = mock_response

        # Call the method
        zenhub_provider, mock_query = mock_zenhub_provider_query
        result = zenhub_provider.get_issues(pipeline_id, repo_ids)

        # Assert the expected response
        assert result == mock_response["data"]["searchIssuesByPipeline"]["nodes"]
        assert '"pipeline_id"' in mock_query.call_args_list[0][0][0]
        assert '["repo_id"]' in mock_query.call_args_list[0][0][0]

    def test_close_issues(self, mock_zenhub_provider_query):
        # Prepare mock data and expected response
        issue_ids = ["issue_id"]
        mock_response = {"successCount": 1}
        _, mock_query = mock_zenhub_provider_query
        mock_query.return_value = mock_response

        # Call the method
        zenhub_provider, mock_query = mock_zenhub_provider_query
        result = zenhub_provider.close_issues(issue_ids)

        # Assert the expected response
        assert result == mock_response
        assert '["issue_id"]' in mock_query.call_args_list[0][0][0]

    def test_query_error_handling(self):
        # Test error handling when _query method raises an exception
        zenhub_provider = ZenHubProvider("fake_token", "fake_url")
        fake_response = Response()
        fake_response.status_code = 500
        fake_session = Mock()
        fake_session.post.return_value = fake_response
        zenhub_provider.session = fake_session

        # Call the method and assert it raises an exception
        with pytest.raises(requests.HTTPError):
            zenhub_provider._query("fake_query")


def test_format_issue():
    issue = {"number": 101, "title": "Sample Issue"}
    assert format_issue(issue) == "#101 (Sample Issue)"


class TestParseWorkspace:
    def test_positive(self):
        node = {"name": "test", "id": "123"}
        workspace = [node]
        assert parse_workspace(workspace, "test") == node

    def test_missing_workspace(self):
        workspace = [{"name": "test", "id": "123"}]
        with pytest.raises(ValueError):
            parse_workspace(workspace, "missing")

    def test_duplicate_workspace(self):
        """The first workspace_name found should be returned if there are duplicates."""
        first_node = {"name": "test", "id": "123"}
        workspace = [first_node, {"name": "test", "id": "456"}]
        assert parse_workspace(workspace, "test") == first_node


class TestParsePipeline:
    def test_positive(self):
        workspace = {"pipelinesConnection": {"nodes": [{"name": "Ready for Prod", "id": "123"}]}}
        assert parse_pipeline(workspace, "Ready for Prod") == "123"

    def test_missing_pipeline(self):
        workspace = {"pipelinesConnection": {"nodes": [{"name": "Ready for Prod", "id": "123"}]}}
        with pytest.raises(ValueError):
            parse_pipeline(workspace, "Missing Pipeline")


class TestParseRepoIds:
    def test_positive(self):
        workspace = {
            "repositoriesConnection": {"nodes": [{"name": "repo1", "id": "123"}, {"name": "repo2", "id": "456"}]}
        }
        repo_names = ["repo1", "repo2"]
        assert parse_repo_ids(workspace, repo_names) == ["123", "456"]

    def test_missing_repo(self):
        workspace = {"repositoriesConnection": {"nodes": [{"name": "repo1", "id": "123"}]}}
        repo_names = ["repo1", "repo2"]
        with pytest.raises(ValueError):
            parse_repo_ids(workspace, repo_names)

    def test_extra_repos(self):
        workspace = {
            "repositoriesConnection": {"nodes": [{"name": "repo1", "id": "123"}, {"name": "repo2", "id": "456"}]}
        }
        repo_names = ["repo1"]
        assert parse_repo_ids(workspace, repo_names) == ["123"]


class TestFilterIssues:
    def test_positive(self, sample_issues):
        issues_to_close, blocked_issues = filter_issues(sample_issues)

        assert len(issues_to_close) == 1
        assert len(blocked_issues) == 0
        assert issues_to_close == [(1, format_issue(sample_issues[0]))]

    def test_with_blocking_issues_in_issue_list(self, sample_issues, blocking_issue):
        sample_issues[0]["blockingIssues"]["nodes"] = [blocking_issue]
        new_issue = sample_issues[0].copy()
        new_issue.update(**blocking_issue)
        sample_issues.append(new_issue)  # add a new issue that is blocking the first issue

        issues_to_close, blocked_issues = filter_issues(sample_issues)

        assert len(issues_to_close) == 2
        assert len(blocked_issues) == 0
        assert issues_to_close == [(1, format_issue(sample_issues[0])), (3, format_issue(blocking_issue))]

    def test_filter_issues_with_blocking_issues(self, sample_issues, blocking_issue):
        sample_issues[0]["blockingIssues"]["nodes"] = [blocking_issue]

        issues_to_close, blocked_issues = filter_issues(sample_issues)

        assert len(issues_to_close) == 0
        assert len(blocked_issues) == 1
        # Validate the structure and content of blocked_issues
        assert blocked_issues[0]["blocking_issues"] == [format_issue(blocking_issue)]
        assert blocked_issues[0]["issue"] == format_issue(sample_issues[0])
        assert blocked_issues[0]["open_prs"] == []

    def test_filter_issues_with_openprs(self, sample_issues, blocking_issue):
        sample_issues[0]["connectedPrs"]["nodes"] = [blocking_issue]

        issues_to_close, blocked_issues = filter_issues(sample_issues)
        assert len(issues_to_close) == 0
        assert len(blocked_issues) == 1
        # Validate the structure and content of blocked_issues
        assert blocked_issues[0]["blocking_issues"] == []
        assert blocked_issues[0]["issue"] == format_issue(sample_issues[0])
        assert blocked_issues[0]["open_prs"] == [format_issue(blocking_issue)]


class TestCloseReadyForProd:
    def test_close_ready_for_prod_no_issues(self, caplog, mock_zenhub_provider):
        caplog.set_level("INFO")
        mock_zenhub_provider.close_issues = Mock(return_value={"data": {"closeIssues": {"successCount": 0}}})

        close_ready_for_prod("example_workspace", ["repo_name"], "pipeline_name", mock_zenhub_provider)

        assert "No issues to close." in caplog.text

    def test_close_ready_for_prod_with_blocked_issues(self, caplog, mock_zenhub_provider):
        caplog.set_level("INFO")
        blocked_issues_data = [
            {
                "issue": "#102 (Sample Issue 2)",
                "open_prs": ["#202 (PR Title 2)"],
                "blocking_issues": ["#301 (Blocked Issue)"],
            }
        ]
        mock_zenhub_provider.get_issues.return_value = [{"id": 1, "title": "Sample Issue"}]
        with patch("scripts.zenhub_close.filter_issues") as mock_filter_issues:
            mock_filter_issues.return_value = ([], blocked_issues_data)

            close_ready_for_prod("example_workspace", ["repo_name"], "pipeline_name", mock_zenhub_provider)

        assert "Not closed: #102 (Sample Issue 2)" in caplog.text
        assert "\tBlocked by:" in caplog.text
        assert "\t\t#301 (Blocked Issue)" in caplog.text
        assert "\tOpen PRs:" in caplog.text
        assert "\t\t#202 (PR Title 2)" in caplog.text

    def test_close_ready_for_prod_with_issues_to_close(self, caplog, mock_zenhub_provider):
        caplog.set_level("INFO")
        issues_to_close_data = [(1, "#101 (Sample Issue 1)")]

        mock_zenhub_provider.get_issues.return_value = [{"id": 1, "title": "Sample Issue"}]
        mock_zenhub_provider.close_issues.return_value = {"data": {"closeIssues": {"successCount": 1}}}
        with patch("scripts.zenhub_close.filter_issues") as mock_filter_issues:
            mock_filter_issues.return_value = (issues_to_close_data, [])

            close_ready_for_prod("example_workspace", ["repo_name"], "pipeline_name", mock_zenhub_provider)

        assert "Closing issues:" in caplog.text
        assert "\t#101 (Sample Issue 1)" in caplog.text
