import logging
import time
from logging import LogRecord
from logging.config import dictConfig
from typing import List

from pythonjsonlogger import jsonlogger
from pythonjsonlogger.jsonlogger import merge_record_extra


# The fields to log using the json logger.
LOGGED_FIELDS = ['levelname', 'asctime', 'name', 'message']
LOG_FORMAT = ' '.join([f"%({field})" for field in LOGGED_FIELDS])


class JsonFormatter(jsonlogger.JsonFormatter):
    default_time_format = '%Y-%m-%dT%H:%M:%S'
    default_msec_format = '%s.%03dZ'

    converter = time.gmtime

    def add_required_fields(self, fields: List[str]) -> None:
        """
        Add additional required fields to to be written in log messages. New fields will be added to the end of the
        `required_fields` list in the order specified by `fields`.

        :param fields: an ordered list of required fields to write to logs.
        :return:
        """

        self._required_fields += [field for field in fields if field not in self._required_fields]
        self._skip_fields = dict(zip(self._required_fields,
                                     self._required_fields))

    def set_required_fields(self, fields: List[str]) -> None:
        """
        Sets the required fields in the order specified in `fields`. Required fields appears in the logs in the order
        listed in `required_fields`.

        :param fields: an ordered list of fields to set `required_fields`
        :return:
        """
        self._required_fields = fields
        self._skip_fields = dict(zip(self._required_fields,
                                     self._required_fields))

    def add_fields(self, log_record: dict, record: LogRecord, message_dict: dict) -> None:
        """
        Adds additional log information from `log_record` to `records. If a required field does not exist in the
        `log_record` then it is not included in the `record`.

        :param log_record: additional fields to add to the `record`.
        :param record: the logRecord to add additional fields too.
        :param message_dict: the log message and extra fields to add to `records`.
        :return:
        """
        for field in self._required_fields:
            value = record.__dict__.get(field)
            if value:
                log_record[field] = value
        log_record.update(message_dict)
        merge_record_extra(record, log_record, reserved=self._skip_fields)


gunicorn_logger = logging.getLogger("gunicorn.error")

dictConfig({
    'version': 1,
    'formatters': {'default': {
        'format': LOG_FORMAT,
        '()': JsonFormatter,
    }},
    'handlers': {
        'wsgi': {
            'class': 'logging.StreamHandler',
            'stream': 'ext://flask.logging.wsgi_errors_stream',
            'formatter': 'default'
        },
        'exceptions':{
            'class': 'logging.StreamHandler',
            'stream': 'ext://sys.stderr',
            'formatter': 'default',
            'level': 'ERROR'
        }
    },
    'root': {
        'level': gunicorn_logger.level,
        'handlers': ['wsgi', 'exceptions']
    }
}
)
