from datetime import datetime, timedelta, timezone
from enum import Enum
from json import JSONEncoder

from backend.common.entities.entity import Entity
from backend.common.corpora_orm import Base


time_zone_info = datetime.now(timezone.utc).astimezone().tzinfo  # Get AWS env time zone info


class CustomJSONEncoder(JSONEncoder):
    "Add support for serializing DateTime, Enums, and SQLAlchemy Base types into JSON"

    def default(self, obj):
        if isinstance(obj, timedelta):
            return str(obj)
        elif isinstance(obj, datetime):
            return obj.timestamp()
        elif isinstance(obj, Enum):
            return str(obj.name)
        elif isinstance(obj, (Base, Entity)):
            return obj.to_dict()
        else:
            return super().default(obj)


class CurationJSONEncoder(CustomJSONEncoder):
    "Add support for serializing DateTime into isoformat"

    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.replace(tzinfo=time_zone_info).isoformat()
        else:
            return super().default(obj)
