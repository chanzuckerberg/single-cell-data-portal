from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.contributor import (
    Contributor,
)

from dcp_prototype.backend.wrangling.migrations.utils.util import merge_dictionary_into
import logging


class Project:
    def __init__(self):
        self.corresponding_old_id = None
        self.project_short_name = None
        self.publication_title = None
        self.publication_doi = None
        self.external_accessions = None

        self.contributors = []

    def populate_from_dcp_one_json_data_frame(self, row):
        self.corresponding_old_id = row.get("provenance.document_id")
        self.project_short_name = row.get("project_core.project_short_name")
        self.publication_title = row.get("publications.0.title")
        self.publication_doi = row.get("publications.0.doi")
        self.external_accessions = row.get("geo_series_accessions.0")

        # Log in case there are multiple publication
        if row.get("publications.1.title"):
            logging.log(
                logging.WARNING,
                f"WARNING: Project {self.project_short_name} contains multiple "
                f"publications but the migration has only captured one of the them "
                f"{self.publication_title}!",
            )
        if row.get("geo_series_accessions.0"):
            logging.log(
                logging.WARNING,
                f"WARNING: Project {self.project_short_name} contains multiple "
                f"external accessions but the migration has only captured one of them "
                f"{self.external_accessions}!",
            )

        contributor_index = 0
        while row.get("contributors." + str(contributor_index) + ".name"):
            contributor = Contributor()
            contributor.populate_from_dcp_one_json_data_frame(row, contributor_index)
            self.add_associated_contributor(contributor)
            contributor_index += 1

    def add_associated_contributor(self, contributor: Contributor):
        if contributor not in self.contributors:
            self.contributors.append(contributor)

    def to_dictionary(self):
        dictionary_representation = {}
        for contributor in self.contributors:
            for key, value in self.__dict__.items():
                if key == "contributors":
                    merge_dictionary_into(
                        dictionary_representation,
                        {key: contributor.name},
                    )
                else:
                    merge_dictionary_into(dictionary_representation, {key: value})
        return dictionary_representation
