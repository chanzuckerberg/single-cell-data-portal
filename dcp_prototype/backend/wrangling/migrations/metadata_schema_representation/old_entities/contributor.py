from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.project import (
    Project,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.project import (
    Project,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.contributor import (
    Contributor,
)
from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_transformer,
)


class Contributor:
    def __init__(self):
        self.name = None
        self.email = None
        self.phone_number = None
        self.corresponding_contributor = None
        self.lab = None
        self.street_address = None
        self.country = None
        self.contributor_role_ontology = None
        self.orcid_id = None

        self.project = None

    def set_project(self, project: Project):
        self.project = project

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

    def convert_to_new_entity(self):
        new_contributor = Contributor()

        new_contributor.name = self.name
        new_contributor.email = self.email
        new_contributor.phone_number = self.email
        new_contributor.corresponding_contributor = self.corresponding_contributor
        new_contributor.lab = self.lab
        new_contributor.street_address = self.street_address
        new_contributor.country = self.country
        new_contributor.contributor_role_ontology = self.contributor_role_ontology
        new_contributor.orcid_id = self.orcid_id
