from typing import List, Optional


class BusinessException(Exception):
    # TODO maybe useful as a base class, maybe not
    pass


# TODO: possibly merge the following 3 classes? They all refer to the same action (validating the collection metadata)
class CollectionCreationException(BusinessException):
    def __init__(self, errors: Optional[List[str]] = None) -> None:
        self.errors: Optional[List[str]] = errors
        super().__init__()


class CollectionUpdateException(BusinessException):
    def __init__(self, errors: Optional[List[str]] = None) -> None:
        self.errors: Optional[List[str]] = errors
        super().__init__()


class CollectionNotFoundException(CollectionUpdateException):
    """
    Raised when a collection is expected to be found, but does not exist
    """

    pass


class CollectionIsPublishedException(CollectionUpdateException):
    """
    Raised when a mutable operation is performed on a published exception
    """

    pass


class InvalidLinkException(BusinessException):
    def __init__(self, errors: Optional[List[str]] = None) -> None:
        self.errors: Optional[List[str]] = errors
        super().__init__()


class CollectionVersionException(BusinessException):
    pass


class ImmutableCollectionException(BusinessException):
    pass


class DatasetIngestException(BusinessException):
    pass


class InvalidURIException(DatasetIngestException):
    """
    Raised when trying to ingest a dataset with an invalid URI
    """

    pass


class MaxFileSizeExceededException(DatasetIngestException):
    """
    Raised when trying to ingest a dataset that is too big
    """

    pass


class DatasetInWrongStatusException(DatasetIngestException):
    """
    Raised when a dataset cannot be updated due to being in a wrong processing status
    """

    pass


class DatasetNotFoundException(BusinessException):
    """
    Raised when a write operation was called on a dataset, but it was not found
    """

    pass


class DatasetUpdateException(BusinessException):
    pass


class CollectionPublishException(BusinessException):
    pass


class ArtifactNotFoundException(BusinessException):
    pass
