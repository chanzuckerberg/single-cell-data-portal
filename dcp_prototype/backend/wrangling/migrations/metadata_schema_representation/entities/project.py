from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_generator,
)

from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.contributor import (
    Contributor,
)
from dcp_prototype.backend.wrangling.migrations.utils.util import merge_dictionary_into


class Project:
    def __init__(self):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.project_short_name = None
        self.publication_title = None
        self.publication_doi = None
        self.external_accessions = None

        self.contributors = []

    def add_associated_contributor(self, contributor: Contributor):
        self.contributors.append(contributor)

    def populate_from_dcp_one_data_row(self, row):
        self.project_short_name = row.get("project.project_core.project_short_name")
        self.publication_title = row.get("project.publications.title")
        self.publication_doi = row.get("project.publications.doi")
        self.external_accessions = row.get("project.geo_series_accessions")

    def to_dictionary(self):
        dictionary_representation = {}
        for contributor in self.contributors:
            for key, value in self.__dict__.items():
                if key == "contributors":
                    merge_dictionary_into(
                        dictionary_representation, {key: contributor.id}
                    )
                else:
                    merge_dictionary_into(dictionary_representation, {key: value})
        return dictionary_representation
