import datetime
from enum import Enum
from json import JSONEncoder


class CustomJSONEncoder(JSONEncoder):
    "Add support for serializing DateTime and Enums types into JSON"

    def default(self, obj):
        if isinstance(obj, datetime.timedelta):
            return str(obj)
        elif isinstance(obj, datetime.datetime):
            return obj.timestamp()
        if isinstance(obj, Enum):
            return str(obj.value)
        else:
            return super().default(self, obj)
