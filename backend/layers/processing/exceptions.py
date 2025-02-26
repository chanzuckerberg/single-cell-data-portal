from dataclasses import dataclass
from typing import List

from backend.layers.common.entities import DatasetStatusKey


class ProcessingException(Exception):
    pass


class ProcessingCanceled(ProcessingException):
    pass


@dataclass
class ValidationAnndataFailed(ProcessingException):
    def __init__(self, errors: List[str]):
        self.errors = errors


class ValidationAtacFailed(ProcessingException):
    def __init__(self, errors: List[str]):
        self.errors = errors


class AddLabelsFailed(ProcessingException):
    failed_status: DatasetStatusKey = DatasetStatusKey.H5AD


class ProcessingFailed(ProcessingException):
    pass


class UploadFailed(ProcessingException):
    pass


@dataclass
class ConversionFailed(ProcessingException):
    failed_status: DatasetStatusKey
