from dataclasses import dataclass
from typing import List

from backend.layers.common.entities import DatasetStatusKey


class ProcessingException(Exception):
    pass


class ProcessingCanceled(ProcessingException):
    pass


@dataclass
class ValidationFailed(ProcessingException):
    errors: List[str]


class ProcessingFailed(ProcessingException):
    pass


class UploadFailed(ProcessingException):
    pass


@dataclass
class ConversionFailed(ProcessingException):
    failed_status: DatasetStatusKey
