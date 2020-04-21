"""
Loads a JSON artifact from S3 or local file
and populates the provided SQL-compatible database.
"""

import json
import os
import sys

import boto3

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from scripts.mock import mock_data

from code.common.browser_orm import (
    DBSessionMaker,
    Project,
    File,
    LibraryConstructionMethod,
    Organ,
    Species,
    Contributor,
    DataRepository,
    ExternalAccession,
    LibraryConstructionMethodJoinProject,
    OrganJoinProject,
    SpeciesJoinProject,
    ContributorJoinProject,
)


def load_from_artifact(session, path_to_file=None):
    if not path_to_file:
        s3 = boto3.client("s3")
        s3.download_file("dcp-test-artifacts", "Artifact.Mar18.json", "artifact.json")
        path_to_file = "artifact.json"

    with open(path_to_file, "r") as f:
        data = json.load(f)
        # augment with mock data
        for project in data["projects"]:
            for key in mock_data[project["id"]]:
                project[key] = mock_data[project["id"]][key]

    organs = {}
    species = {}
    libraries = {}
    contributors = {}

    # populate data repository
    data_repo_names = ["Array Express", "INSDC Project", "GEO Series", "Biostudies"]
    session.add_all([DataRepository(name=name) for name in data_repo_names])
    session.commit()

    # populate project
    for project in data["projects"]:
        session.add(
            Project(
                id=project["id"],
                title=project["title"],
                label=project["label"],
                description=project["description"],
                biosample_categories=",".join(project["biosample_categories"]),
                development_stages=",".join(project["donor_development_stages_at_collection"]),
                diseases=",".join(project["donor_diseases"]),
                cell_isolation_methods=",".join(project["cell_isolation_methods"]),
                cell_types=",".join(project["selected_cell_types"]),
                cell_count=project["cell_count"],
                paired_end=",".join([str(e) for e in project["paired_end"]]),
                nucleic_acid_sources=",".join(project["nucleic_acid_sources"]),
                input_nucleic_acid_molecules=",".join(project["input_nucleic_acid_molecules"]),
                publication_title=project["publication_title"],
                publication_doi=project["publication_doi"],
                cxg_enabled=project["cxg_enabled"],
            )
        )
        session.commit()

        # populate organ + project join
        for o in project["organs"]:
            if o not in organs:
                organs[o] = len(organs) + 1
                session.add(Organ(name=o))
            session.add(OrganJoinProject(organ_id=organs[o], project_id=project["id"]))

        # populate species + project join
        for s in project["donor_species"]:
            if s not in species:
                species[s] = len(species) + 1
                session.add(Species(name=s))
            session.add(SpeciesJoinProject(species_id=species[s], project_id=project["id"]))

        # populate library_construction_method + project join
        for l in project["library_construction_methods"]:
            if l not in libraries:
                libraries[l] = len(libraries) + 1
                session.add(LibraryConstructionMethod(name=l))
            session.add(
                LibraryConstructionMethodJoinProject(
                    library_construction_method_id=libraries[l], project_id=project["id"]
                )
            )

        # populate contributors + project join
        for c in project["contributors"]:
            if c["name"] not in contributors:
                contributors[c["name"]] = {"key": len(contributors) + 1, "institution": c["institution"]}

                names = c["name"].split(",")
                if len(names) == 3:
                    first = names[0]
                    middle = names[1]
                    last = names[2]
                elif len(names) == 2:
                    first = names[0]
                    middle = ""
                    last = names[1]
                else:
                    first, last = c["name"].split(" ")
                    middle = ""

                session.add(
                    Contributor(first_name=first, middle_name=middle, last_name=last, institution=c["institution"])
                )
            session.add(ContributorJoinProject(contributor_id=contributors[c["name"]]["key"], project_id=project["id"]))

        # populate external accessions
        for accession in project["array_express_accessions"]:
            session.add(ExternalAccession(project_id=project["id"], data_repository_id=1, accession=accession))
        for accession in project["insdc_project_accessions"]:
            session.add(ExternalAccession(project_id=project["id"], data_repository_id=2, accession=accession))
        for accession in project["geo_series_accessions"]:
            session.add(ExternalAccession(project_id=project["id"], data_repository_id=3, accession=accession))
        for accession in project["biostudies_accessions"]:
            session.add(ExternalAccession(project_id=project["id"], data_repository_id=4, accession=accession))
        session.commit()

    # file
    session.add_all(
        File(
            id=file["id"],
            project_id=file["project_id"],
            filename=file["filename"],
            file_format=file["file_format"],
            file_size=file["file_size"],
            file_type=file["type"],
            s3_uri=file["s3_uri"],
        )
        for file in data["files"]
    )

    session.commit()
    session.close()


if __name__ == "__main__":
    load_from_artifact(DBSessionMaker().session())
