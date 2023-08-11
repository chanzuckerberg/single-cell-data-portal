import requests  # type: ignore
from connexion.exceptions import ProblemException


class ServerErrorHTTPException(ProblemException):
    def __init__(
        self, detail: str = "An internal server error has occurred. Please try again later.", *args, **kwargs
    ) -> None:
        super().__init__(
            *args,
            **kwargs,
            status=requests.codes.server_error,
            title="Internal Server Error",
            detail=detail,
        )


class UnauthorizedError(ProblemException):
    def __init__(self, detail: str = "Invalid Authentication Credentials", *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, status=401, title="Invalid Authentication Credentials", detail=detail)


class ExpiredCredentialsError(ProblemException):
    def __init__(self, detail: str = "Invalid Authentication Credentials", *args, **kwargs) -> None:
        super().__init__(
            *args,
            **kwargs,
            status=401,
            title="Invalid Authentication Credentials",
            detail="Your authentication credentials have expired. Please provide updated "
            "credentials before trying again.",
        )


class TooLargeHTTPException(ProblemException):
    def __init__(self, detail: str = "The file requested is too large to process.", *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, status=413, title="Request Entity Too Large", detail=detail)


class InvalidParametersHTTPException(ProblemException):
    _default_detail = "One or more parameters is invalid."

    def __init__(self, detail: str = None, *args, **kwargs) -> None:  # type: ignore
        detail = detail if detail else self._default_detail
        super().__init__(*args, **kwargs, status=400, title="Bad Request", detail=detail)


class ForbiddenHTTPException(ProblemException):
    def __init__(self, detail: str = "User is not authorized to access this resource.", *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, status=requests.codes.forbidden, title="Forbidden", detail=detail)


class NotFoundHTTPException(ProblemException):
    def __init__(self, detail: str = "Resource not found.", *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, status=requests.codes.not_found, title="Not Found", detail=detail)


class MethodNotAllowedException(ProblemException):
    def __init__(self, detail: str = "Method not allowed.", *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, status=requests.codes.not_allowed, title="Not Allowed", detail=detail)


class ConflictException(ProblemException):
    def __init__(self, detail: str = "A duplicate resource already exists.", *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, status=requests.codes.conflict, title="Conflict", detail=detail)


class GoneHTTPException(ProblemException):
    def __init__(self, detail: str = "Resource has been removed", *args, **kwargs) -> None:
        super().__init__(*args, **kwargs, status=requests.codes.gone, title="Gone", detail=detail)
