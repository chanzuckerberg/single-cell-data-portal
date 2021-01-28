import typing
from datetime import datetime

from sqlalchemy import and_

from ..utils.db_utils import clone
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

        # Setting Defaults
        links = links if links else []

        new_db_object = DbCollection(
            visibility=visibility,
            name=name,
            description=description,
            owner=owner,
            contact_name=contact_name,
            contact_email=contact_email,
            data_submission_policy_version=data_submission_policy_version,
            **kwargs,
        )

        new_db_object.links = [
            DbCollectionLink(collection_id=new_db_object.id, collection_visibility=visibility, **link) for link in links
        ]
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
        Return a collection if the user is the owner of a collection.

        :param collection_uuid: the uuid of the collection
        :param visibility: the visibility of the collection
        :param user: the uuid of the user.
        :return: a collection if the user is the owner of the collection else None
        """
        filters = [cls.table.id == collection_uuid, cls.table.owner == user, cls.table.visibility == visibility]
        collection = cls.db.session.query(cls.table).filter(*filters).one_or_none()
        return cls(collection) if collection else None

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

    _hide = ["user", "owner", "tombstone"]

    def reshape_for_api(self) -> dict:
        """
        Reshape the collection to match the expected api output.
        :return: A dictionary that can be converted into JSON matching the expected api response.
        """

        result = self.to_dict(remove_none=True)
        # Reshape the data to match.
        for hidden in self._hide:
            result.pop(hidden, None)
        result["links"] = [
            dict(link_url=link["link_url"], link_name=link["link_name"] or "", link_type=link["link_type"])
            for link in result["links"]
        ]

        result["datasets"] = [ds for ds in result["datasets"] if ds.get("tombtone") is not True]
        for dataset in result["datasets"]:
            dataset["dataset_deployments"] = dataset.pop("deployment_directories")
            dataset["dataset_assets"] = dataset.pop("artifacts")
            dataset.pop("tombstone")

        return result

    def publish(self):
        """
        Given a private collection, set the collection to public.

        """
        # Create a public collection with the same uuid and same fields
        public_collection = clone(self.db_object, primary_key=dict(id=self.id, visibility=CollectionVisibility.PUBLIC))
        self.db.session.add(public_collection)
        # Copy over relationships
        for link in self.links:
            link.collection_visibility = CollectionVisibility.PUBLIC
        for dataset in self.datasets:
            dataset.collection_visibility = CollectionVisibility.PUBLIC
        self.db.session.commit()
        self.delete()
        self.db_object = public_collection
