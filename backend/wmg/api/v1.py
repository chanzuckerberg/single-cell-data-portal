from uuid import UUID, uuid4

from flask import jsonify


def latest_snapshot() -> UUID:
    return jsonify(snapshot_id=uuid4().hex)
