import typing
import uuid

from .entity import Entity
from ..corpora_orm import DbDataset, DbDatasetArtifact, DbDeploymentDirectory, DbContributor, DbDatasetContributor


class Dataset(Entity):
    table = DbDataset

    def __init__(self, db_object: DbDataset):
        super().__init__(db_object)

    @classmethod
    def create(
        cls,
        revision: int = 0,
        name: str = "",
        organism: str = "",
        organism_ontology: str = "",
        tissue: str = "",
        tissue_ontology: str = "",
        assay: str = "",
        assay_ontology: str = "",
        disease: str = "",
        disease_ontology: str = "",
        sex: str = "",
        ethnicity: str = "",
        ethnicity_ontology: str = "",
        source_data_location: str = "",
        preprint_doi: str = "",
        publication_doi: str = "",
        artifacts: list = None,
        contributors: list = None,
        deployment_directories: list = None,
        **kwargs,
    ) -> "Dataset":
        """
        Creates a new dataset and related objects and store in the database. UUIDs are generated for all new table
        entries.
        """
        primary_key = str(uuid.uuid4())

        # Setting Defaults
        artifacts = artifacts if artifacts else []
        deployment_directories = deployment_directories if deployment_directories else []
        contributors = contributors if contributors else []

        new_db_object = DbDataset(
            id=primary_key,
            revision=revision,
            name=name,
            organism=organism,
            organism_ontology=organism_ontology,
            tissue=tissue,
            tissue_ontology=tissue_ontology,
            assay=assay,
            assay_ontology=assay_ontology,
            disease=disease,
            disease_ontology=disease_ontology,
            sex=sex,
            ethnicity=ethnicity,
            ethnicity_ontology=ethnicity_ontology,
            source_data_location=source_data_location,
            preprint_doi=preprint_doi,
            publication_doi=publication_doi,
            artifacts=cls._create_sub_objects(artifacts, DbDatasetArtifact, add_columns=dict(dataset_id=primary_key)),
            deployment_directories=cls._create_sub_objects(
                deployment_directories, DbDeploymentDirectory, add_columns=dict(dataset_id=primary_key),
            ),
            **kwargs,
        )

        #  Linking many contributors to many datasets
        add_contributors = cls.create_contributors(contributors)
        contributor_dataset_ids = [
            dict(contributor_id=contributor.id, dataset_id=primary_key) for contributor in add_contributors
        ]
        dataset_contributor = cls._create_sub_objects(contributor_dataset_ids, DbDatasetContributor)

        cls.db.session.add(new_db_object)
        cls.db.session.flush()
        cls.db.session.add_all(dataset_contributor)
        cls.db.commit()

        return cls(new_db_object)

    @classmethod
    def create_contributors(cls, contributors: typing.List[dict]) -> typing.List[DbContributor]:
        """
        Take a list of contributor. Retrieve them from the database if they exists, otherwise create them. Return a list
        of DbContributors.
        :param contributors: a dictionary containing all of the columns for a DbContributor row.
        :return: a list of DbContributor objects.
        """
        existing_contributors = []
        new_contributors = []
        for contributor in contributors:
            email = contributor.get("email")
            if email:
                result = cls.db.query([DbContributor], [DbContributor.email == email])
                if result:
                    existing_contributors.append(result[0])
                else:
                    new_contributors.append(contributor)

        #  Create new contributors
        new_contributors = cls._create_sub_objects(new_contributors, DbContributor)
        if new_contributors:
            cls.db.session.add_all(new_contributors)
            cls.db.commit()
        return [*new_contributors, *existing_contributors]
