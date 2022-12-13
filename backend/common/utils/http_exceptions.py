import requests
from connexion.exceptions import ProblemException


class ServerErrorHTTPException(ProblemException):
    def __init__(
        self, detail: str = "An internal server error has occurred. Please try again later.", *args, **kwargs
    ) -> None:
        super().__init__(
            status=requests.codes.server_error, title="Internal Server Error", detail=detail, *args, **kwargs
        )


class UnauthorizedError(ProblemException):
    def __init__(self, detail: str = "Invalid Authentication Credentials", *args, **kwargs) -> None:
        super().__init__(status=401, title="Invalid Authentication Credentials", detail=detail, *args, **kwargs)


class ExpiredCredentialsError(ProblemException):
    def __init__(self, detail: str = "Invalid Authentication Credentials", *args, **kwargs) -> None:
        super().__init__(
            status=401,
            title="Invalid Authentication Credentials",
            detail="Your authentication credentials have expired. Please provide updated "
            "credentials before trying again.",
            *args,
            **kwargs,
        )


class TooLargeHTTPException(ProblemException):
    def __init__(self, detail: str = "The file requested is too large to process.", *args, **kwargs) -> None:
        super().__init__(status=413, title="Request Entity Too Large", detail=detail, *args, **kwargs)


class InvalidParametersHTTPException(ProblemException):
    _default_detail = "One or more parameters is invalid."

    def __init__(self, detail: str = None, *args, **kwargs) -> None:
        detail = detail if detail else self._default_detail
        super().__init__(status=400, title="Bad Request", detail=detail, *args, **kwargs)


class ForbiddenHTTPException(ProblemException):
    def __init__(self, detail: str = "User is not authorized to access this resource.", *args, **kwargs) -> None:
        super().__init__(status=requests.codes.forbidden, title="Forbidden", detail=detail, *args, **kwargs)


class NotFoundHTTPException(ProblemException):
    def __init__(self, detail: str = "Resource not found.", *args, **kwargs) -> None:
        super().__init__(status=requests.codes.not_found, title="Not Found", detail=detail, *args, **kwargs)


class MethodNotAllowedException(ProblemException):
    def __init__(self, detail: str = "Method not allowed.", *args, **kwargs) -> None:
        super().__init__(status=requests.codes.not_allowed, title="Not Allowed", detail=detail, *args, **kwargs)


class ConflictException(ProblemException):
    def __init__(self, detail: str = "A duplicate resource already exists.", *args, **kwargs) -> None:
        super().__init__(status=requests.codes.conflict, title="Conflict", detail=detail, *args, **kwargs)


class GoneHTTPException(ProblemException):
    def __init__(self, detail: str = "Resource has been removed", *args, **kwargs) -> None:
        super().__init__(status=requests.codes.gone, title="Gone", detail=detail, *args, **kwargs)
