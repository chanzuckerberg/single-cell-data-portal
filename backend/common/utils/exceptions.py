class CorporaException(Exception):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class MaxFileSizeExceededException(CorporaException):
    pass


class InvalidFileFormatException(CorporaException):
    pass


class NonExistentCollectionException(CorporaException):
    pass


class InvalidProcessingStateException(CorporaException):
    pass


class NonExistentDatasetException(CorporaException):
    pass


class MarkerGeneCalculationException(Exception):
    pass
