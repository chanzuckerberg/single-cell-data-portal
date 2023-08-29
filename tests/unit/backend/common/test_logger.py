import pytest

from backend.common.logging_config import LogSuppressed


def test_log_suppressed__caught(caplog):
    # check a log message is written when an exception is suppressed
    with LogSuppressed(ValueError):
        raise ValueError("Test Exception")
    assert "Suppressed Exception" in caplog.text


def test_log_suppressed__missed():
    # the exception passes through
    with pytest.raises(ValueError), LogSuppressed(KeyError):
        raise ValueError("Test Exception")
