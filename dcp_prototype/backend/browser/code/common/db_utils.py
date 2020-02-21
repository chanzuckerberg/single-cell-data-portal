import os
import sys

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.code.common.browser_orm import (
    DBSessionMaker,
    LibraryPrepProtocol,
    Tissue,
    Species,
    LibraryPrepProtocolJoinProject,
    TissueJoinProject,
    SpeciesJoinProject,
)


def _get_project_assays(project_id, session=None):
    """
    Query the DB to return all assays that are represented in a given project.
    :param project_id: Project to return assays for
    :param session: SQLAlchemy DBSession
    :return: list of assay ontology IDs
    """
    if not session:
        session = DBSessionMaker().session()

    assays = []
    for result in (
        session.query(LibraryPrepProtocolJoinProject, LibraryPrepProtocol)
        .filter(
            (LibraryPrepProtocolJoinProject.library_prep_protocol_id == LibraryPrepProtocol.id),
            LibraryPrepProtocolJoinProject.project_id == project_id,
        )
        .all()
    ):
        assays.append(result.LibraryPrepProtocol.construction_method_ontology)

    return assays


def _get_project_tissues(project_id, session=None):
    """
    Query the DB to return all tissues that are represented in a given project.
    :param project_id: Project to return tissues for
    :param session: SQLAlchemy DBSession
    :return: list of tissue ontology IDs
    """
    if not session:
        session = DBSessionMaker().session()

    tissues = []
    for result in (
        session.query(TissueJoinProject, Tissue)
        .filter(TissueJoinProject.tissue_id == Tissue.id, TissueJoinProject.project_id == project_id,)
        .all()
    ):
        tissues.append(result.Tissue.tissue_ontology)

    return tissues


def _get_project_species(project_id, session=None):
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
        species.append(result.Species.species_ontology)

    return species
