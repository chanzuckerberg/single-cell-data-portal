from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_generator,
)

from pandas import Series


class Contributor:
    def __init__(self):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.name = None
        self.email = None
        self.phone_number = None
        self.corresponding_contributor = None
        self.lab = None
        self.street_address = None
        self.country = None
        self.contributor_role_ontology = None
        self.orcid_id = None

    def populate_from_dcp_one_data_row(self, row):
        self.name = row.get("project.contributors.name")
        self.email = row.get("project.contributors.email")
        self.phone_number = row.get("project.contributors.phone")
        corresponding_contributor = row.get(
            "project.contributors.corresponding_contributor"
        )
        self.corresponding_contributor = (
            True
            if corresponding_contributor and corresponding_contributor == "yes"
            else False
        )
        self.lab = row.get("project.contributors.laboratory")
        self.street_address = row.get("project.contributors.address")
        self.country = row.get("project.contributors.country")
        self.contributor_role_ontology = row.get(
            "project.contributors.project_role.ontology_label"
        )
        self.orcid_id = row.get("project.contributors.orcid_id")

    def to_dictionary(self):
        return self.__dict__
