import os
import sys

import chalice
from chalice import Chalice

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "chalicelib"))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from browser.code.rds.browser_orm import (DBSessionMaker,
                                          Project, File,
                                          LibraryPrepProtocol, Tissue, Species,
                                          LibraryPrepProtocolJoinProject, TissueJoinProject,
                                          SpeciesJoinProject)

app = Chalice(app_name='browser-api')
session = DBSessionMaker().session()


@app.route('/')
def index():
    return "render documentation here"


@app.route('/projects')
def get_projects():
    projects = []
    for project in session.query(Project).all():
        projects.append({
            'id': project.id,
            'title': project.title,
            'assays': _get_project_assays(project.id),
            'tissues': _get_project_tissues(project.id),
            'species': _get_project_species(project.id),
            'cell_count': project.cell_count,
        })

    return chalice.Response(status_code=200,
                            headers={
                                'Content-Type': "application/json",
                                'Access-Control-Allow-Origin': "*"
                            },
                            body=projects)


@app.route('/projects/{project_id}')
def get_project(project_id):
    project = session.query(Project).get(project_id)
    response_body = {
        'id': project.id,
        'title': project.title,
        'label': project.label,
        'assays': _get_project_assays(project.id),
        'tissues': _get_project_tissues(project.id),
        'species': _get_project_species(project.id),
        'description': project.description,
        'category': project.category,
        'developmental_stage': project.developmental_stage,
        'disease_ontology': project.disease_ontology,
        'sample_type': project.sample_type,
        'organ_part': project.organ_part,
        'analysis_protocol': project.analysis_protocol.split(','),
        'cell_count': project.cell_count,
        'donor_count': project.donor_count,
        'publication_title': project.publication_title,
        'publication_doi': project.publication_doi,
        'contact_name': project.contact_name,
        'contact_institution': project.contact_institution,
        'contact_email': project.contact_email
    }

    return chalice.Response(status_code=200,
                            headers={
                                'Content-Type': "application/json",
                                'Access-Control-Allow-Origin': "*"
                            },
                            body=response_body)


@app.route('/projects/{project_id}/files')
def get_project_files(project_id):
    files = []
    for file in session.query(File).filter(File.project_id == project_id).limit(100):
        files.append({
            'id': file.id,
            'filename': file.filename,
            'file_format': file.file_format,
            'file_size': file.file_size,
            'species': file.species,
            'library_construction_method_ontology':
                file.library_construction_method_ontology,
            'tissue_ontology': file.tissue_ontology
        })

    return chalice.Response(status_code=200,
                            headers={
                                'Content-Type': "application/json",
                                'Access-Control-Allow-Origin': "*"
                            },
                            body=files)


@app.route('/files/{file_id}')
def get_file(file_id):
    file = session.query(File).get(file_id)
    response_body = {
        's3_uri': file.s3_uri
    }
    return chalice.Response(status_code=200,
                            headers={
                                'Content-Type': "application/json",
                                'Access-Control-Allow-Origin': "*"
                            },
                            body=response_body)


def _get_project_assays(project_id):
    assays = []
    for result in session.query(
        LibraryPrepProtocolJoinProject, LibraryPrepProtocol
    ).filter(
        (LibraryPrepProtocolJoinProject.library_prep_protocol_id ==
         LibraryPrepProtocol.id),
        LibraryPrepProtocolJoinProject.project_id == project_id
    ).all():
        assays.append(result.LibraryPrepProtocol.construction_method_ontology)

    return assays


def _get_project_tissues(project_id):
    tissues = []
    for result in session.query(TissueJoinProject, Tissue).filter(
            TissueJoinProject.tissue_id == Tissue.id,
            TissueJoinProject.project_id == project_id).all():
        tissues.append(result.Tissue.tissue_ontology)

    return tissues


def _get_project_species(project_id):
    species = []
    for result in session.query(SpeciesJoinProject, Species).filter(
            SpeciesJoinProject.species_id == Species.id,
            SpeciesJoinProject.project_id == project_id).all():
        species.append(result.Species.species_ontology)

    return species
