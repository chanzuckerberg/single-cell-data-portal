class CorporaException(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class AuthorizationError(Exception):
    status_code = 403
    msg = "The user is not authorized to access this resource"

    def __init__(self, message: str = None, status_code: int = None, payload: dict = None):
        Exception.__init__(self)
        self.message = message if message else self.msg
        if status_code is not None:
            self.status_code = status_code
        self.payload = payload

    def to_dict(self):
        rv = dict(self.payload or ())
        rv["message"] = self.message
        return rv
