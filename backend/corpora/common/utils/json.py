import datetime
from enum import Enum
from json import JSONEncoder

from ..entities.entity import Entity
from ..corpora_orm import Base


class CustomJSONEncoder(JSONEncoder):
    "Add support for serializing DateTime and Enums types into JSON"

    def default(self, obj):
        if isinstance(obj, datetime.timedelta):
            return str(obj)
        elif isinstance(obj, datetime.datetime):
            return obj.timestamp()
        elif isinstance(obj, Enum):
            return str(obj.name)
        elif isinstance(obj, (Base, Entity)):
            return obj.to_dict()
        else:
            return super().default(self, obj)
