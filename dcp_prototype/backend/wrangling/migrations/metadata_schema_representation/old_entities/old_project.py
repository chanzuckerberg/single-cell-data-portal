from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_contributor import (
    OldContributor,
)
import logging
from copy import deepcopy

from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_transformer,
)
from dcp_prototype.backend.ledger.code.common.ledger_orm import Project


class OldProject:
    def __init__(self):
        self.corresponding_old_id = None
        self.project_short_name = None
        self.publication_title = None
        self.publication_doi = None
        self.external_accessions = []

    def populate_from_dcp_one_json_data_frame(self, row):
        self.corresponding_old_id = row.get("provenance.document_id")
        self.project_short_name = row.get("project_core.project_short_name")
        self.publication_title = row.get("publications.0.title")
        self.publication_doi = row.get("publications.0.doi")

        accession_index = 0
        while row.get(f"geo_series_accessions.{str(accession_index)}"):
            self.external_accessions.append(
                row.get(f"geo_series_accessions.{str(accession_index)}")
            )
        accession_index = 0
        while row.get(f"insdc_project_accessions.{str(accession_index)}"):
            self.external_accessions.append(
                row.get(f"insdc_project_accessions.{str(accession_index)}")
            )

        # Log in case there are multiple publications
        if row.get("publications.1.title"):
            logging.log(
                logging.WARNING,
                f"WARNING: Project {self.project_short_name} contains multiple "
                f"publications but the migration has only captured one of the them "
                f"{self.publication_title}!",
            )

        contributors = []
        contributor_index = 0
        while row.get("contributors." + str(contributor_index) + ".name"):
            contributor = OldContributor()
            contributor.populate_from_dcp_one_json_data_frame(row, contributor_index)
            contributors.append(contributor)
            contributor_index += 1

        return contributors

    def to_dictionary(self):
        return deepcopy(self.__dict__)

    def convert_to_new_entity(self):
        project_id = hca_accession_transformer(
            Project.__name__, self.corresponding_old_id
        )

        project = Project(
            id=project_id,
            project_short_name=self.project_short_name,
            publication_title=self.publication_title,
            publication_doi=self.publication_doi,
            external_accessions=self.external_accessions,
        )

        return project
