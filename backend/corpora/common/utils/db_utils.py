import os
import typing

from ..corpora_orm import (
    Base,
    DBSessionMaker,
    Project,
    File,
    LibraryConstructionMethod,
    Organ,
    Species,
    Contributor,
    ContributorJoinProject,
    LibraryConstructionMethodJoinProject,
    OrganJoinProject,
    SpeciesJoinProject,
)


class DbUtils:
    def __init__(self):
        self.session = DBSessionMaker().session()
        self.engine = self.session.get_bind()

    def _get(self, table: Base, entity_id: str) -> typing.Union[Base, None]:
        """
        Query a table row by its primary key
        :param table: SQLAlchemy Table to query
        :param entity_id: Primary key of desired row
        :return: SQLAlchemy Table object, None if not found
        """
        return self.session.query(table).get(entity_id)

    def _query(self, table_args: typing.List[Base], filter_args: typing.List[bool] = None) -> typing.List[Base]:
        """
        Query the database using the current DB session
        :param table_args: List of SQLAlchemy Tables to query/join
        :param filter_args: List of SQLAlchemy filter conditions
        :return: List of SQLAlchemy query response objects
        """
        return (
            self.session.query(*table_args).filter(*filter_args).all()
            if filter_args
            else self.session.query(*table_args).all()
        )

    @staticmethod
    def _parse_multivalue(value: str) -> typing.List[str]:
        """
        Parses a CSV string representing multiple values into a list
        :param value: Comma-separated value to parse
        :return: List of strings
        """
        return value.split(",") if value else []

    def _is_test_db(self) -> bool:
        """
        Tests whether the current DB connection
        is to a test database or not
        :return: True if test DB, else False
        """
        return self.engine.driver == "pysqlite"

    def create(self) -> None:
        """
        Drop and recreate all tables.
        This operation is only supported for the test SQLite database.
        Use the admin tool to perform this action on deployed environments.
        """
        if not self._is_test_db():
            raise EnvironmentError(
                f"{os.environ['DEPLOYMENT_STAGE']} is not a test environment. "
                f"Operation not supported. Please use the admin tool to perform this action."
            )

        Base.metadata.drop_all(bind=self.engine)
        Base.metadata.create_all(bind=self.engine)

    def query_projects(self) -> typing.List[dict]:
        """
        Query the DB for all projects
        :return: List of project metadata dicts
        """
        projects = [
            {
                "id": project.id,
                "title": project.title,
                "assays": self.query_project_assays(project.id),
                "organs": self.query_project_organs(project.id),
                "species": self.query_project_species(project.id),
                "cell_count": project.cell_count,
            }
            for project in self._query([Project])
        ]

        return projects

    def query_project(self, project_id: str) -> typing.Union[dict, None]:
        """
        Query the DB for a project by its project ID
        :param project_id: Project ID
        :return: Project metadata dict
        """
        project = self._get(Project, project_id)

        return (
            {
                "id": project.id,
                "title": project.title,
                "label": project.label,
                "assays": self.query_project_assays(project.id),
                "organs": self.query_project_organs(project.id),
                "species": self.query_project_species(project.id),
                "contributors": self.query_project_contributors(project.id),
                "description": project.description,
                "biosample_categories": self._parse_multivalue(project.biosample_categories),
                "development_stages": self._parse_multivalue(project.development_stages),
                "diseases": self._parse_multivalue(project.diseases),
                "cell_isolation_methods": self._parse_multivalue(project.cell_isolation_methods),
                "cell_types": self._parse_multivalue(project.cell_types),
                "cell_count": project.cell_count,
                "paired_end": self._parse_multivalue(project.paired_end),
                "nucleic_acid_sources": self._parse_multivalue(project.nucleic_acid_sources),
                "input_nucleic_acid_molecules": self._parse_multivalue(project.input_nucleic_acid_molecules),
                "publication_title": project.publication_title,
                "publication_doi": project.publication_doi,
                "cxg_enabled": project.cxg_enabled,
            }
            if project
            else None
        )

    def query_file(self, file_id: str) -> typing.Union[dict, None]:
        """
        Query the DB for a file by its file ID
        :param file_id: File ID
        :return: File query result
        """
        file = self._get(File, file_id)

        return (
            {
                "id": file.id,
                "project_id": file.project_id,
                "filename": file.filename,
                "file_format": file.file_format,
                "file_size": file.file_size,
                "file_type": file.file_type,
                "s3_uri": file.s3_uri,
            }
            if file
            else None
        )

    def query_project_assays(self, project_id: str) -> typing.List[str]:
        """
        Query the DB to return all assays that are represented in a given project.
        :param project_id: Project ID to return assays for
        :return: list of assay names
        """
        assays = []
        for result in self._query(
            table_args=[LibraryConstructionMethodJoinProject, LibraryConstructionMethod],
            filter_args=[
                LibraryConstructionMethodJoinProject.library_construction_method_id == LibraryConstructionMethod.id,
                LibraryConstructionMethodJoinProject.project_id == project_id,
            ],
        ):
            assays.append(result.LibraryConstructionMethod.name)

        return assays

    def query_project_organs(self, project_id: str) -> typing.List[str]:
        """
        Query the DB to return all organs that are represented in a given project.
        :param project_id: Project ID to return organs for
        :return: list of organ names
        """
        organs = []
        for result in self._query(
            table_args=[OrganJoinProject, Organ],
            filter_args=[OrganJoinProject.organ_id == Organ.id, OrganJoinProject.project_id == project_id],
        ):
            organs.append(result.Organ.name)

        return organs

    def query_project_species(self, project_id: str) -> typing.List[str]:
        """
        Query the DB to return all species that are represented in a given project.
        :param project_id: Project ID to return species for
        :return: list of species labels
        """
        species = []
        for result in self._query(
            table_args=[SpeciesJoinProject, Species],
            filter_args=[SpeciesJoinProject.species_id == Species.id, SpeciesJoinProject.project_id == project_id],
        ):
            species.append(result.Species.name)

        return species

    def query_project_contributors(self, project_id: str) -> typing.List[dict]:
        """
        Query the DB to return all contributors associated with a given project.
        :param project_id: Project ID to return contributors for
        :return: list of contributors
        """
        contributor_query = self._query(
            table_args=[Contributor, ContributorJoinProject],
            filter_args=[
                Contributor.id == ContributorJoinProject.contributor_id,
                ContributorJoinProject.project_id == project_id,
            ],
        )

        contributors = [
            {
                "name": f"{result.Contributor.first_name} {result.Contributor.last_name}",
                "institution": f"{result.Contributor.institution}",
            }
            for result in contributor_query
        ]

        return contributors

    def query_downloadable_project_files(self, project_id: str) -> typing.List[dict]:
        """
        Query the DB to return all downloadable files for a project.
        :param project_id: Project to return files for
        :return: list of file metadata objects
        """
        files = []
        for file in self._query(
            table_args=[File], filter_args=[File.project_id == project_id, File.file_format == "LOOM"]
        ):
            files.append(
                {
                    "id": file.id,
                    "filename": file.filename,
                    "file_format": file.file_format,
                    "file_type": file.file_type,
                    "file_size": file.file_size,
                }
            )

        return files
