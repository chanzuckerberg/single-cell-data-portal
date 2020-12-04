import requests
from connexion.exceptions import ProblemException


class CorporaException(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class ServerErrorHTTPException(ProblemException):
    def __init__(
        self, detail: str = "An internal server error has occurred. Please try again later.", *args, **kwargs
    ) -> None:
        super().__init__(
            status=requests.codes.server_error, title="Internal Server Error", detail=detail, *args, **kwargs
        )


class TooLargeHTTPException(ProblemException):
    def __init__(self, detail: str = "The file requested is too large to process.", *args, **kwargs) -> None:
        super().__init__(status=413, title="Request Entity Too Large", detail=detail, *args, **kwargs)


class InvalidParametersHTTPException(ProblemException):
    def __init__(self, detail: str = "One or more parameters is invalid.", *args, **kwargs) -> None:
        super().__init__(status=400, title="Bad Request", detail=detail, *args, **kwargs)


class ForbiddenHTTPException(ProblemException):
    def __init__(self, detail: str = "User is not authorized to access this resource.", *args, **kwargs) -> None:
        super().__init__(status=requests.codes.forbidden, title="Forbidden", detail=detail, *args, **kwargs)


class NotFoundHTTPException(ProblemException):
    def __init__(self, detail: str = "Resource not found.", *args, **kwargs) -> None:
        super().__init__(status=requests.codes.not_found, title="Not Found", detail=detail, *args, **kwargs)
