from dataclasses import dataclass
from typing import List

from backend.layers.common.entities import DatasetStatusKey


class ProcessingException(Exception):
    pass


class ProcessingCanceled(ProcessingException):
    pass


@dataclass
class ValidationAnndataFailed(ProcessingException):
    errors: List[str]


class ValidationAtacFailed(ProcessingException):
    errors: List[str]


class AddLabelsFailed(ProcessingException):
    failed_status: DatasetStatusKey = DatasetStatusKey.H5AD


class ProcessingFailed(ProcessingException):
    pass


class UploadFailed(ProcessingException):
    pass


@dataclass
class ConversionFailed(ProcessingException):
    failed_status: DatasetStatusKey
