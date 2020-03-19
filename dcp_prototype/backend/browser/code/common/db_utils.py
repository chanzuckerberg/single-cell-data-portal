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


def get_project_assays(project_id, session=None):
    """
    Query the DB to return all assays that are represented in a given project.
    :param project_id: Project to return assays for
    :param session: SQLAlchemy DBSession
    :return: list of assay names
    """
    if not session:
        session = DBSessionMaker().session()

    assays = []
    for result in (
        session.query(LibraryConstructionMethodJoinProject, LibraryConstructionMethod)
        .filter(
            (LibraryConstructionMethodJoinProject.library_construction_method_id == LibraryConstructionMethod.id),
            LibraryConstructionMethodJoinProject.project_id == project_id
        )
        .all()
    ):
        assays.append(result.LibraryConstructionMethod.name)

    return assays


def get_project_organs(project_id, session=None):
    """
    Query the DB to return all organs that are represented in a given project.
    :param project_id: Project to return organs for
    :param session: SQLAlchemy DBSession
    :return: list of organ names
    """
    if not session:
        session = DBSessionMaker().session()

    organs = []
    for result in (
        session.query(OrganJoinProject, Organ)
        .filter(OrganJoinProject.organ_id == Organ.id, OrganJoinProject.project_id == project_id)
        .all()
    ):
        organs.append(result.Organ.name)

    return organs


def get_project_species(project_id, session=None):
    """
    Query the DB to return all species that are represented in a given project.
    :param project_id: Project to return species for
    :param session: SQLAlchemy DBSession
    :return: list of species labels
    """
    if not session:
        session = DBSessionMaker().session()

    species = []
    for result in (
        session.query(SpeciesJoinProject, Species)
        .filter(SpeciesJoinProject.species_id == Species.id, SpeciesJoinProject.project_id == project_id,)
        .all()
    ):
        species.append(result.Species.name)

    return species


def get_project_contributors(project_id, session=None):
    if not session:
        session = DBSessionMaker().session()

    contributor_query = session.query(Contributor, ContributorJoinProject).filter(
        Contributor.id == ContributorJoinProject.contributor_id,
        ContributorJoinProject.project_id == project_id
    ).all()

    contributors = [{
        'name': f"{result.Contributor.first_name} {result.Contributor.last_name}",
        'institution': f"{result.Contributor.institution}"
    } for result in contributor_query]

    return contributors


def get_downloadable_project_files(project_id, session=None):
    """
    Query the DB to return all downloadable files for a project.
    :param project_id: Project to return files for
    :param session: SQLAlchemy DBSession
    :return: list of file metadata objects
    """
    if not session:
        session = DBSessionMaker().session()

    files = []
    for file in session.query(File).filter(File.project_id == project_id).filter(File.file_format == "LOOM"):
        files.append(
            {
                "id": file.id,
                "filename": file.filename,
                "file_format": file.file_format,
                "file_type": file.file_type,
                "file_size": file.file_size
            }
        )

    return files
