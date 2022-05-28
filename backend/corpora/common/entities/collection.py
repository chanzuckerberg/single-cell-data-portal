import typing
from datetime import datetime
from sqlalchemy import and_
from sqlalchemy.orm import Session

from . import Dataset
from .entity import Entity
from .geneset import Geneset
from .collection_link import CollectionLink
from ..corpora_orm import (
    CollectionLinkType,
    DbCollection,
    DbCollectionLink,
    CollectionVisibility,
    generate_uuid,
)
from ..utils.db_helpers import clone


class Collection(Entity):
    table = DbCollection
    list_attributes = (DbCollection.id, DbCollection.created_at)

    def __init__(self, db_object: DbCollection):
        super().__init__(db_object)

    @classmethod
    def create(
        cls,
        session: Session,
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
        new_db_object.links = [DbCollectionLink(collection_id=new_db_object.id, **link) for link in links]
        session.add(new_db_object)
        session.commit()
        return cls(new_db_object)

    @classmethod
    def get_collection(
        cls,
        session: Session,
        collection_uuid: str = None,
        visibility: CollectionVisibility = None,
        revision_of: str = None,
        include_tombstones: bool = False,
        owner: typing.Optional[str] = None,
    ) -> typing.Union["Collection", None]:
        """
        Given the collection_uuid, retrieve a live collection.
        :param session: the database session object.
        :param collection_uuid:
        :param visibility: the visibility of the collection
        :param revision_of: the id of the associated public collection iff this collection is a revision
        :param include_tombstones: If true, the collection is returned even if it has been tombstoned.
        :param owner: A user id use to check if the user is the owner of the collection. If the user id matches the
        owner then the collection is returned. If this parameters is not included then owner is not used as a filter.
        :return: the collection if it matches the filter.
        """
        filters = []
        if visibility:
            filters.append(cls.table.visibility == visibility)
        if collection_uuid:
            filters.append(cls.table.id == collection_uuid)
        if revision_of:
            filters.append(cls.table.revision_of == revision_of)
        if owner:
            filters.append(cls.table.owner == owner)
        if not include_tombstones:
            filters.append(cls.table.tombstone == False)  # noqa
        collection = session.query(cls.table).filter(*filters).one_or_none()

        return cls(collection) if collection else None

    @classmethod
    def list_collections_in_time_range(cls, session: Session, *args, **kwargs):
        return cls.list_attributes_in_time_range(
            session, *args, filters=[DbCollection.visibility == CollectionVisibility.PUBLIC.name], **kwargs
        )

    @classmethod
    def list_attributes_in_time_range(
        cls,
        session: Session,
        to_date: int = None,
        from_date: int = None,
        filters: list = None,
        list_attributes: list = None,
    ) -> typing.List[typing.Dict]:
        """
        Queries the database for Entities that have been created within the specified time range. Return only the
        entity attributes in `list_attributes`.

        :param session: The SQLAlchemy database Session
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

    @classmethod
    def list_public_datasets_for_index(cls, session: Session) -> typing.List[typing.Dict]:
        """
        Return a list of all the datasets and associated metadata. For efficiency reasons, this only returns the fields
        inside the `dataset` table and doesn't include relationships.
        """

        attrs = [
            Dataset.table.id,
            Dataset.table.name,
            Dataset.table.collection_id,
            Dataset.table.tissue,
            Dataset.table.disease,
            Dataset.table.assay,
            Dataset.table.organism,
            Dataset.table.cell_count,
            Dataset.table.cell_type,
            Dataset.table.sex,
            Dataset.table.ethnicity,
            Dataset.table.development_stage,
            Dataset.table.is_primary_data,
            Dataset.table.mean_genes_per_cell,
            Dataset.table.schema_version,  # Required for schema manipulation
            Dataset.table.explorer_url,
            Dataset.table.published_at,
            Dataset.table.revised_at,
        ]

        def to_dict(db_object):
            _result = {}
            for _field in db_object._fields:
                _value = getattr(db_object, _field)
                if _value is None:
                    continue
                _result[_field] = getattr(db_object, _field)
            return _result

        filters = [~Dataset.table.tombstone, cls.table.visibility == CollectionVisibility.PUBLIC.name]

        results = [
            to_dict(result)
            for result in session.query(Dataset.table).join(cls.table).filter(*filters).with_entities(*attrs).all()
        ]

        for result in results:
            Dataset.transform_organism_for_schema_2_0_0(result)
            Dataset.transform_sex_for_schema_2_0_0(result)
            Dataset.enrich_development_stage_with_ancestors(result)

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
        is_existing_collection = False
        if self.revision_of:
            # This is a revision of a published Collection
            public_collection = Collection.get_collection(self.session, collection_uuid=self.revision_of)
            revision = self.to_dict(
                remove_attr=("updated_at", "created_at", "visibility", "id", "revision_of", "published_at"),
                remove_relationships=True,
            )
            revision["data_submission_policy_version"] = data_submission_policy_version
            public_collection.update(  # This function call deletes old links (keep_links param defaults to False)
                commit=False,
                **revision,
            )
            for link in self.links:
                CollectionLink(link).update(collection_id=self.revision_of, commit=False)
            is_existing_collection = True

        else:
            # This is the first publishing of this Collection
            self.update(
                commit=False,
                visibility=CollectionVisibility.PUBLIC,
                published_at=now,
                data_submission_policy_version=data_submission_policy_version,
                keep_links=True,
            )

        has_dataset_changes = False
        for dataset in self.datasets:
            revised_dataset = Dataset(dataset)
            original = Dataset.get(self.session, revised_dataset.original_id) if revised_dataset.original_id else None
            if original and public_collection.check_has_dataset(original):
                dataset_is_changed = original.publish_revision(revised_dataset, now)
                if dataset_is_changed:
                    has_dataset_changes = True
            else:
                # The dataset is new
                revised_dataset.publish_new(now)
                has_dataset_changes = True

        self.session.flush()
        self.session.expire_all()

        if is_existing_collection:
            if has_dataset_changes:
                public_collection.update(commit=False, remove_attr="revised_at", revised_at=now, keep_links=True)
            self.delete()
            self.db_object = public_collection.db_object

        self.session.commit()

    def create_revision(self) -> "Collection":
        """
        Generate a collection revision from a public collection
        :return: collection revision.

        """
        revision_collection = clone(
            self.db_object,
            primary_key=dict(id=generate_uuid()),
            visibility=CollectionVisibility.PRIVATE,
            revision_of=self.id,
        )
        self.session.add(revision_collection)
        for link in self.links:
            self.session.add(clone(link, collection_id=revision_collection.id))
        for dataset in self.datasets:
            Dataset(dataset).create_revision(revision_collection.id)
        self.session.commit()
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

    def update(self, links: list = None, keep_links=False, **kwargs) -> None:
        """
        Update an existing collection to match provided the parameters. The specified columns are replaced.
        :param links: links to create and connect to the collection. If present, the existing attached entries will
         be removed and replaced with new entries.
        :param keep_links: boolean - whether or not links need to be preserved. Links are preserved if True.
        :param kwargs: Any other fields in the dataset that will be replaced.
        """
        links = links if links else []
        if not keep_links:
            for link in self.links:
                self.session.delete(link)

        new_objs = [DbCollectionLink(collection_id=self.id, **link) for link in links]
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
        return self.id == dataset.collection_id

    def get_doi(self) -> str:
        doi = [link for link in self.links if link.link_type == CollectionLinkType.DOI]
        if doi:
            return doi[0].link_url
        else:
            return None
