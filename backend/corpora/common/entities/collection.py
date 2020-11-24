import typing
import uuid
from datetime import datetime

from sqlalchemy import and_

from .entity import Entity
from ..corpora_orm import DbCollection, DbCollectionLink, CollectionVisibility


class Collection(Entity):
    table = DbCollection
    list_attributes = (DbCollection.id, DbCollection.created_at)

    def __init__(self, db_object: DbCollection):
        super().__init__(db_object)

    @classmethod
    def create(
        cls,
        visibility: CollectionVisibility,
        name: str = "",
        description: str = "",
        owner: str = "",
        contact_name: str = "",
        contact_email: str = "",
        links: list = None,
        data_submission_policy_version: str = "",
        **kwargs,
    ) -> "Collection":
        """
        Create a new Collection and related objects and store in the database. UUIDs are generated for all new table
        entries.
        """
        primary_key = str(uuid.uuid4())

        # Setting Defaults
        links = links if links else []

        new_db_object = DbCollection(
            id=primary_key,
            visibility=visibility,
            name=name,
            description=description,
            owner=owner,
            contact_name=contact_name,
            contact_email=contact_email,
            data_submission_policy_version=data_submission_policy_version,
            links=cls._create_sub_objects(
                links, DbCollectionLink, add_columns=dict(collection_id=primary_key, collection_visibility=visibility)
            ),
            **kwargs,
        )

        cls.db.session.add(new_db_object)
        cls.db.commit()
        return cls(new_db_object)

    @classmethod
    def get_collection(cls, collection_uuid, visibility=CollectionVisibility.PUBLIC.name):
        """
        Given the collection_uuid, retrieve a live collection.
        :param collection_uuid:
        """
        return cls.get((collection_uuid, visibility))

    @classmethod
    def if_owner(
        cls, collection_uuid: str, visibility: CollectionVisibility, user: str
    ) -> typing.Union[DbCollection, None]:
        """

        :param collection_uuid: the uuid of the collection
        :param visibility: the visibility of the collection
        :param user: the uuid of the user.
        :return: a collection if the user is the owner of the collection else None
        """
        filters = [cls.table.id == collection_uuid, cls.table.owner == user, cls.table.visibility == visibility]
        return cls.db.session.query(cls.table).filter(*filters).one_or_none()

    @classmethod
    def list_collections_in_time_range(cls, *args, **kwargs):
        return cls.list_attributes_in_time_range(
            *args, filters=[DbCollection.visibility == CollectionVisibility.PUBLIC.name], **kwargs
        )

    @classmethod
    def list_attributes_in_time_range(
        cls, to_date: int = None, from_date: int = None, filters: list = None, list_attributes: list = None
    ) -> typing.List[typing.Dict]:
        """
        Queries the database for Entities that have been created within the specified time range. Return only the
        entity attributes in `list_attributes`.

        :param to_date: If provided, only lists collections that were created before this date. Format of param is Unix
        timestamp since the epoch in UTC timezone.
        :param from_date: If provided, only lists collections that were created after this date. Format of param is Unix
        timestamp since the epoch in UTC timezone.
        :param filters: additional filters to apply to the query.
        :param list_attributes: A list of entity attributes to return. If None, the class default is used.
        :return: The results is a list of flattened dictionaries containing the `list_attributes`
        """

        filters = filters if filters else []
        list_attributes = list_attributes if list_attributes else cls.list_attributes
        table = cls.table

        def to_dict(db_object):
            _result = {}
            for _field in db_object._fields:
                _result[_field] = getattr(db_object, _field)
            return _result

        if to_date:
            filters.append(cls.table.created_at <= datetime.fromtimestamp(to_date))
        if from_date:
            filters.append(table.created_at >= datetime.fromtimestamp(from_date))

        results = [
            to_dict(result)
            for result in cls.db.session.query(table).with_entities(*list_attributes).filter(and_(*filters)).all()
        ]

        return results

    @classmethod
    def get_submission(cls, collection_uuid):
        """
        Given the collection_uuid, retrieve a live collection.
        :param collection_uuid:
        """
        return cls.get((collection_uuid, CollectionVisibility.PRIVATE.name))

    @classmethod
    def list_submissions(cls, *args, **kwargs):
        return cls.list_attributes_in_time_range(
            *args,
            filters=[DbCollection.visibility == CollectionVisibility.PRIVATE.name],
            list_attributes=[
                DbCollection.id,
                DbCollection.name,
                DbCollection.owner,
            ],
            **kwargs,
        )

    def reshape_for_api(self) -> dict:
        """
        Reshape the collection to match the expected api output.
        :return: A dictionary that can be converted into JSON matching the expected api response.
        """
        result = self.to_dict()
        # Reshape the data to match.
        result.pop("user", None)
        result.pop("owner", None)
        result["links"] = [
            dict(link_url=link["link_url"], link_name=link["link_name"] or "", link_type=link["link_type"])
            for link in result["links"]
        ]
        for dataset in result["datasets"]:
            dataset["dataset_deployments"] = dataset.pop("deployment_directories")
            dataset["dataset_assets"] = dataset.pop("artifacts")

        return result
