import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.code.common.browser_orm import (
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

    def _get(self, table, entity_id):  # noqa
        return self.session.query(table).get(entity_id)

    def _query(self, table_args, filter_args=None):  # noqa
        return (
            self.session.query(*table_args).filter(*filter_args).all()
            if filter_args
            else self.session.query(*table_args).all()
        )

    def query_projects(self):
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

    def query_project(self, project_id):
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
                "biosample_categories": project.biosample_categories.split(","),
                "development_stages": project.development_stages.split(","),
                "diseases": project.diseases.split(","),
                "cell_isolation_methods": project.cell_isolation_methods.split(","),
                "cell_types": project.cell_types.split(","),
                "cell_count": project.cell_count,
                "paired_end": project.paired_end.split(","),
                "nucleic_acid_sources": project.nucleic_acid_sources.split(","),
                "input_nucleic_acid_molecules": project.input_nucleic_acid_molecules.split(","),
                "publication_title": project.publication_title,
                "publication_doi": project.publication_doi,
                "cxg_enabled": project.cxg_enabled,
            }
            if project
            else None
        )

    def query_file(self, file_id):
        """
        Query the DB for a file by its file ID
        :param file_id: File ID
        :return: File query result
        """
        return self._get(File, file_id)

    def query_project_assays(self, project_id):
        """
        Query the DB to return all assays that are represented in a given project.
        :param project_id: Project ID to return assays for
        :return: list of assay names
        """
        assays = []
        for result in self._query(
            table_args=[LibraryConstructionMethodJoinProject, LibraryConstructionMethod],
            filter_args=[
                (LibraryConstructionMethodJoinProject.library_construction_method_id == LibraryConstructionMethod.id),
                LibraryConstructionMethodJoinProject.project_id == project_id,
            ],
        ):
            assays.append(result.LibraryConstructionMethod.name)

        return assays

    def query_project_organs(self, project_id):
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

    def query_project_species(self, project_id):
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

    def query_project_contributors(self, project_id):
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

    def query_downloadable_project_files(self, project_id):
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
