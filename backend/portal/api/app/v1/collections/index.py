from flask import make_response, jsonify, g

from backend.api_server.db import dbconnect
from backend.common.corpora_orm import DbCollection, CollectionVisibility
from backend.common.entities import Collection


@dbconnect
def get():
    # TODO (ebezzi): this is very similar to `get` above. Eventually they should be consolidated
    db_session = g.db_session

    filtered_collection = Collection.list_attributes_in_time_range(
        db_session,
        filters=[DbCollection.visibility == CollectionVisibility.PUBLIC],
        list_attributes=[
            DbCollection.id,
            DbCollection.name,
            DbCollection.published_at,
            DbCollection.revised_at,
            DbCollection.publisher_metadata,
        ],
    )

    # Remove entries where the value is None
    updated_collection = []
    for d in filtered_collection:
        updated_collection.append({k: v for k, v in d.items() if v is not None})

    return make_response(jsonify(updated_collection), 200)
