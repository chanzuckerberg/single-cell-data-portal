from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_generator,
)


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
            if corresponding_contributor
            and (
                corresponding_contributor.uppe() == "YES"
                or corresponding_contributor.upper() is "TRUE"
            )
            else False
        )
        self.lab = row.get("project.contributors.laboratory")
        self.street_address = row.get("project.contributors.address")
        self.country = row.get("project.contributors.country")
        self.contributor_role_ontology = row.get(
            "project.contributors.project_role.ontology_label"
        )
        self.orcid_id = row.get("project.contributors.orcid_id")

    def populate_from_dcp_one_json_data_frame(self, row, index_in_project_row):
        prefix = f"contributors.{str(index_in_project_row)}."
        self.name = row.get(prefix + "name")
        self.email = row.get(prefix + "email")
        self.phone_number = row.get(prefix + "phone")
        corresponding_contributor = row.get(prefix + "corresponding_contributor")
        self.corresponding_contributor = (
            True
            if corresponding_contributor
            and (
                corresponding_contributor.uppe() == "YES"
                or corresponding_contributor.upper() is "TRUE"
            )
            else False
        )
        self.lab = row.get(prefix + "laboratory")
        self.street_address = row.get(prefix + "address")
        self.country = row.get(prefix + "country")
        self.contributor_role_ontology = row.get(prefix + "project_role.ontology")
        self.orcid_id = row.get(prefix + "orcid_id")

    def to_dictionary(self):
        return self.__dict__
