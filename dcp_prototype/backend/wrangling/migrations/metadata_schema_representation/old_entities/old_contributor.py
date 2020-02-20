from dcp_prototype.backend.ledger.code.common.ledger_orm import Contributor
from dcp_prototype.backend.wrangling.migrations.utils.id_generator import hca_accession_generator
from copy import deepcopy


class OldContributor:
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

    def populate_from_dcp_one_json_data_frame(self, row, index_in_project_row):
        prefix = f"contributors.{str(index_in_project_row)}."
        self.name = row.get(prefix + "name")
        self.email = row.get(prefix + "email")
        self.phone_number = row.get(prefix + "phone")
        corresponding_contributor = row.get(prefix + "corresponding_contributor")
        self.corresponding_contributor = (
            True
            if corresponding_contributor
            and (corresponding_contributor.upper() == "YES" or corresponding_contributor.upper() == "TRUE")
            else False
        )
        self.lab = row.get(prefix + "laboratory")
        self.street_address = row.get(prefix + "address")
        self.country = row.get(prefix + "country")
        self.contributor_role_ontology = row.get(prefix + "project_role.ontology")
        self.orcid_id = row.get(prefix + "orcid_id")

    def to_dictionary(self):
        return deepcopy(self.__dict__)

    def convert_to_new_entity(self):
        contributor_id = hca_accession_generator(Contributor.__name__)

        contributor = Contributor(
            id=contributor_id,
            name=self.name,
            email=self.email,
            phone_number=self.phone_number,
            corresponding_contributor=self.corresponding_contributor,
            lab=self.lab,
            street_address=self.street_address,
            country=self.country,
            contributor_role_ontology=self.contributor_role_ontology,
            orcid_id=self.orcid_id,
        )

        return contributor
