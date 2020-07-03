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
    ) -> "Dataset":
        """
        Creates a new dataset and related objects
        :param revision:
        :param name:
        :param organism:
        :param organism_ontology:
        :param tissue:
        :param tissue_ontology:
        :param assay:
        :param assay_ontology:
        :param disease:
        :param disease_ontology:
        :param sex:
        :param ethnicity:
        :param ethnicity_ontology:
        :param source_data_location:
        :param preprint_doi:
        :param publication_doi:
        :param artifacts:
        :param contributors:
        :param deployment_directories:
        :return:
        """
        uuid = cls.db.generate_id(DbDataset)

        """
        Prevent accidentally linking an existing row to a different Dataset. This maintains the relationship of one
        to many for artifacts and deployment_directories
        """
        if artifacts:
            [artifact.pop("id", None) for artifact in artifacts]  # sanitize of ids
        if deployment_directories:
            [deployment_directory.pop("id", None) for deployment_directory in deployment_directories]  # sanitize of ids

        new_db_object = DbDataset(
            id=uuid,
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
            artifacts=cls._create_sub_objects(artifacts, DbDatasetArtifact, dataset_id=uuid),
            deployment_directories=cls._create_sub_objects(
                deployment_directories, DbDeploymentDirectory, dataset_id=uuid
            ),
        )

        #  Linking many contributors to many datasets
        contributors = cls._create_sub_objects(contributors, DbContributor)
        contributor_dataset_ids = [dict(contributor_id=contributor.id, dataset_id=uuid) for contributor in contributors]
        dataset_contributor = cls._create_sub_objects(contributor_dataset_ids, DbDatasetContributor)

        cls.db.session.add(new_db_object)
        cls.db.session.add_all(contributors)
        cls.db.session.commit()

        cls.db.session.add_all(dataset_contributor)
        cls.db.session.commit()

        return cls(new_db_object)
