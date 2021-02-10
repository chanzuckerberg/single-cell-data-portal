from backend.corpora.common.utils.exceptions import CorporaException


class ValidationFailed(Exception):
    pass


class ProcessingFailed(Exception):
    pass


class CorporaTombstoneException(CorporaException):
    pass
