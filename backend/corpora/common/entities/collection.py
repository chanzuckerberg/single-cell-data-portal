import typing
from datetime import datetime
from sqlalchemy import and_

from . import Dataset
from .entity import Entity
from .geneset import Geneset
from ..corpora_orm import CollectionLinkType, DbCollection, DbCollectionLink, CollectionVisibility
from ..utils.db_helpers import clone


class Collection(Entity):
    table = DbCollection
    list_attributes = (DbCollection.id, DbCollection.created_at)

    def __init__(self, db_object: DbCollection):
        super().__init__(db_object)

    @classmethod
    def create(
        cls,
        session,
        visibility: CollectionVisibility,
        name: str = "",
        description: str = "",
        owner: str = "",
        contact_name: str = "",
        contact_email: str = "",
        links: list = None,
        data_submission_policy_version: str = "",
        publisher_metadata: dict = None,
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
            publisher_metadata=publisher_metadata,
            **kwargs,
        )

        new_db_object.links = [
            DbCollectionLink(collection_id=new_db_object.id, collection_visibility=visibility, **link) for link in links
        ]
        session.add(new_db_object)
        session.commit()
        return cls(new_db_object)

    @classmethod
    def get_collection(
        cls,
        session,
        collection_uuid: str,
        visibility: str = CollectionVisibility.PUBLIC.name,
        include_tombstones: bool = False,
        owner: typing.Optional[str] = None,
    ) -> typing.Union["Collection", None]:
        """
        Given the collection_uuid, retrieve a live collection.
        :param session: the database session object.
        :param collection_uuid:
        :param visibility: the visibility of the collection
        :param include_tombstones: If true, the collection is returned even if it has been tombstoned.
        :param owner: A user id use to check if the user is the owner of the collection. If the user id matches the
        owner then the collection is returned. If this parameters is not included then owner is not used as a filter.
        :return: the collection if it matches the filter.
        """
        filters = [cls.table.id == collection_uuid, cls.table.visibility == visibility]
        if owner:
            filters.append(cls.table.owner == owner)
        if not include_tombstones:
            filters.append(cls.table.tombstone == False)  # noqa
        collection = session.query(cls.table).filter(*filters).one_or_none()
        return cls(collection) if collection else None

    @classmethod
    def list_collections_in_time_range(cls, session, *args, **kwargs):
        return cls.list_attributes_in_time_range(
            session, *args, filters=[DbCollection.visibility == CollectionVisibility.PUBLIC.name], **kwargs
        )

    @classmethod
    def list_attributes_in_time_range(
        cls, session, to_date: int = None, from_date: int = None, filters: list = None, list_attributes: list = None
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
        filters.append(cls.table.tombstone == False)  # noqa
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
            for result in session.query(table).with_entities(*list_attributes).filter(and_(*filters)).all()
        ]

        return results

    def reshape_for_api(self, tombstoned_datasets=False) -> dict:
        """
        Reshape the collection to match the expected api output.
        :tombstoned_datasets: Determines if tombtoned datasets will be included.
        :return: A dictionary that can be converted into JSON matching the expected api response.
        """

        result = self.to_dict(remove_none=True)
        # Reshape the data to match.
        for hidden in ["user", "owner", "tombstone"]:
            result.pop(hidden, None)
        result["links"] = [
            dict(link_url=link["link_url"], link_name=link.get("link_name", ""), link_type=link["link_type"])
            for link in result["links"]
        ]
        datasets = []
        for dataset in result["datasets"]:
            dataset["dataset_deployments"] = []
            explorer_url = dataset.pop("explorer_url", None)
            if explorer_url:
                dataset["dataset_deployments"].append({"url": explorer_url})
            dataset["dataset_assets"] = dataset.pop("artifacts")
            if dataset.get("tombstone"):
                if tombstoned_datasets and not dataset["published"]:
                    datasets.append(dataset)
                else:
                    continue
            else:
                datasets.append(dataset)

            Dataset.transform_sex_for_schema_2_0_0(dataset)
            Dataset.transform_organism_for_schema_2_0_0(dataset)

        result["datasets"] = datasets
        return result

    def publish(self, data_submission_policy_version):
        """
        Given a private collection, set the collection to public.
        """
        # Timestamp for published_at and revised_at
        now = datetime.utcnow()
        # Create a public collection with the same uuid and same fields
        public_collection = Collection.get_collection(self.session, self.id, CollectionVisibility.PUBLIC)
        is_existing_collection = False

        if public_collection:
            revision = self.to_dict(
                remove_attr=("updated_at", "created_at", "visibility", "id"), remove_relationships=True
            )
            revision["data_submission_policy_version"] = data_submission_policy_version
            public_collection.update(
                commit=False,
                **revision,
            )
            is_existing_collection = True
        # A published collection with the same uuid does not already exist.
        # This is a new collection.
        else:
            public_collection = Collection(
                clone(
                    self.db_object,
                    primary_key=dict(id=self.id, visibility=CollectionVisibility.PUBLIC),
                    # We want to update published_at only when the collection is first published.
                    published_at=now,
                    data_submission_policy_version=data_submission_policy_version,
                )
            )
            self.session.add(public_collection)

        # Copy over relationships
        for link in self.links:
            link.collection_visibility = CollectionVisibility.PUBLIC

        has_dataset_changes = False
        for dataset in self.datasets:
            revision = Dataset(dataset)
            original = Dataset.get(self.session, revision.original_id) if revision.original_id else None
            if original and public_collection.check_has_dataset(original):
                dataset_is_changed = original.publish_revision(revision, now)
                if dataset_is_changed:
                    has_dataset_changes = True
            else:
                # The dataset is new
                revision.publish_new(now)
                has_dataset_changes = True

        if is_existing_collection and has_dataset_changes:
            public_collection.update(commit=False, remove_attr="revised_at", revised_at=now)

        self.session.commit()
        # commit expires the session, need to retrieve the original private collection to delete it
        private_collection = Collection.get_collection(self.session, self.id, CollectionVisibility.PRIVATE)
        private_collection.delete()
        self.db_object = public_collection.db_object

    def revision(self) -> "Collection":
        """
        Generate a collection revision from a public collection
        :return: collection revision.

        """
        revision_collection = clone(
            self.db_object, primary_key=dict(id=self.id, visibility=CollectionVisibility.PRIVATE)
        )
        self.session.add(revision_collection)
        for link in self.links:
            self.session.add(clone(link, collection_id=self.id, collection_visibility=CollectionVisibility.PRIVATE))
        self.session.commit()
        for dataset in self.datasets:
            Dataset(dataset).create_revision()
        return Collection(revision_collection)

    def tombstone_collection(self):
        self.update(tombstone=True, commit=False)
        for geneset in self.genesets:
            Geneset.get(self.session, geneset.id).delete(commit=False)
        for dataset in self.datasets:
            ds = Dataset.get(self.session, dataset.id, include_tombstones=True)
            ds.asset_deletion()
            ds.tombstone_dataset_and_delete_child_objects()
        self.session.commit()

    def update(self, links: list = None, **kwargs) -> None:
        """
        Update an existing collection to match provided the parameters. The specified columns are replaced.
        :param links: links to create and connect to the collection. If present, the existing attached entries will
         be removed and replaced with new entries.
        :param kwargs: Any other fields in the dataset that will be replaced.
        """
        links = links if links else []
        for link in self.links:
            self.session.delete(link)
        new_objs = [
            DbCollectionLink(collection_id=self.id, collection_visibility=self.visibility, **link) for link in links
        ]
        self.session.add_all(new_objs)

        super().update(**kwargs)

    def delete(self):
        for dataset in self.datasets:
            ds = Dataset(dataset)
            if not ds.published:
                ds.asset_deletion()
            ds.delete()
        super().delete()

    def check_has_dataset(self, dataset: Dataset) -> bool:
        """Check that a dataset is part of the collection"""
        return all([self.id == dataset.collection_id, self.visibility == dataset.collection_visibility])

    def get_doi(self) -> str:
        doi = [link for link in self.links if link.link_type == CollectionLinkType.DOI]
        if doi:
            return doi[0].link_url
        else:
            return None
