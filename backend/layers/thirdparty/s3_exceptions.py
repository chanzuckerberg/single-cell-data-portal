from typing import List, Optional


class S3Exception(Exception):
    # TODO maybe useful as a base class, maybe not
    pass


class S3DeleteException(S3Exception):
    def __init__(self, errors: Optional[List[dict]] = None) -> None:
        self.errors: Optional[List[dict]] = errors
        super().__init__()


class IllegalS3RecursiveDelete(S3Exception):
    pass
