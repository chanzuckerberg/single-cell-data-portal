from unittest.mock import patch

import pytest

from scripts.zenhub_close import close_ready_for_prod, filter_issues, format_issue


@pytest.fixture
def sample_issues():
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
def blocking_issue():
    return {"id": 3, "state": "OPEN", "number": 301, "title": "Blocked Issue"}


def test_format_issue():
    issue = {"number": 101, "title": "Sample Issue"}
    assert format_issue(issue) == "#101 (Sample Issue)"


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
    def test_close_ready_for_prod_no_issues(self, caplog):
        caplog.set_level("INFO")
        with patch("scripts.zenhub_close.get_issues") as mock_get_issues, patch(
            "scripts.zenhub_close.filter_issues"
        ) as mock_filter_issues, patch("scripts.zenhub_close.close_issues"):
            mock_get_issues.return_value = []
            mock_filter_issues.return_value = ([], [])

            close_ready_for_prod()

            assert "No issues to close." in caplog.text

    def test_close_ready_for_prod_with_blocked_issues(self, caplog):
        caplog.set_level("INFO")
        blocked_issues_data = [
            {
                "issue": "#102 (Sample Issue 2)",
                "open_prs": ["#202 (PR Title 2)"],
                "blocking_issues": ["#301 (Blocked Issue)"],
            }
        ]

        with patch("scripts.zenhub_close.get_issues") as mock_get_issues, patch(
            "scripts.zenhub_close.filter_issues"
        ) as mock_filter_issues, patch("scripts.zenhub_close.close_issues"):
            mock_get_issues.return_value = [{"id": 1, "title": "Sample Issue"}]
            mock_filter_issues.return_value = ([], blocked_issues_data)

            close_ready_for_prod()

            assert "Not closed: #102 (Sample Issue 2)" in caplog.text
            assert "\tBlocked by:" in caplog.text
            assert "\t\t#301 (Blocked Issue)" in caplog.text
            assert "\tOpen PRs:" in caplog.text
            assert "\t\t#202 (PR Title 2)" in caplog.text

    def test_close_ready_for_prod_with_issues_to_close(self, caplog):
        caplog.set_level("INFO")
        issues_to_close_data = [(1, "#101 (Sample Issue 1)")]

        with patch("scripts.zenhub_close.get_issues") as mock_get_issues, patch(
            "scripts.zenhub_close.filter_issues"
        ) as mock_filter_issues, patch("scripts.zenhub_close.close_issues") as mock_close_issues:
            mock_get_issues.return_value = [{"id": 1, "title": "Sample Issue"}]
            mock_filter_issues.return_value = (issues_to_close_data, [])
            mock_close_issues.return_value = {"data": {"closeIssues": {"successCount": 1}}}

            close_ready_for_prod()

            assert "Closing issues:" in caplog.text
            assert "\t#101 (Sample Issue 1)" in caplog.text
