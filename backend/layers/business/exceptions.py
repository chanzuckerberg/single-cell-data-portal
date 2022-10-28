class BusinessException(Exception):
    # TODO maybe useful as a base class, maybe not
    pass


class CollectionCreationException(BusinessException):
    pass


class CollectionUpdateException(BusinessException):
    pass


class InvalidLinkException(BusinessException):
    pass


class DatasetIngestException(BusinessException):
    pass


class CollectionPublishException(BusinessException):
    pass
