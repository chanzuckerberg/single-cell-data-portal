from typing import List, Optional


class BusinessException(Exception):
    # TODO maybe useful as a base class, maybe not
    pass


# TODO: possibly merge the following 3 classes? They all refer to the same action (validating the collection metadata)
class CollectionCreationException(BusinessException):
    
    errors: Optional[List[str]]

    def __init__(self, errors: Optional[List[str]] = None) -> None:
        self.errors = errors
        super().__init__()


class CollectionUpdateException(BusinessException):

    errors: Optional[List[str]]

    def __init__(self, errors: Optional[List[str]] = None) -> None:
        self.errors = errors
        super().__init__()


class InvalidLinkException(BusinessException):
    
    errors: Optional[List[str]]

    def __init__(self, errors: Optional[List[str]] = None) -> None:
        self.errors = errors
        super().__init__()

class ImmutableCollectionException(BusinessException):
    pass


class DatasetIngestException(BusinessException):
    pass

class DatasetUpdateException(BusinessException):
    pass


class CollectionPublishException(BusinessException):
    pass

class ArtifactNotFoundException(BusinessException):
    pass