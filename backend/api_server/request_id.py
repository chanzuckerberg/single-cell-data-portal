"""
Inspired by: http://blog.mcpolemic.com/2016/01/18/adding-request-ids-to-flask.html
"""
import logging
import uuid

import flask


def generate_request_id():
    return uuid.uuid4()


def get_request_id():
    return getattr(flask.g, "request_id", None)


class RequestIdFilter(logging.Filter):
    # This is a logging filter that makes the request ID available for use in
    # the logging format. Note that we're checking if we're in a request
    # context, as we may want to log things before Flask is fully loaded.
    def filter(self, record):
        record.request_id = flask.g.request_id if flask.has_request_context() else ""
        return True
