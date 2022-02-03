from flask import make_response


def dummy():
    return make_response(dict(message="I found your gene!"))
