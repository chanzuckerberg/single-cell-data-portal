import requests
from connexion.exceptions import ProblemException


class CorporaException(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class ForbiddenHTTPException(ProblemException):
    def __init__(self, detail: str = "User is not authorized to access this resource.", *args, **kwargs) -> None:
        super().__init__(status=requests.codes.forbidden, title="Forbidden", detail=detail, *args, **kwargs)
