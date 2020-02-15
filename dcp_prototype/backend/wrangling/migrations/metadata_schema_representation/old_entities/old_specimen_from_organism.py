from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_donor_organism import (  # noqa
    OldDonorOrganism,
)

from copy import deepcopy
import logging


class OldSpecimenFromOrganism:
    def __init__(self):
        self.corresponding_old_id = None
        self.organ = None

        self.donor_organism = None

    def populate_from_dcp_one_json_data_frame(self, row):
        self.corresponding_old_id = row.get("provenance.document_id")
        self.organ = row.get("organ.ontology")

    def set_donor_organism(self, donor_organism: OldDonorOrganism):
        if (
            self.donor_organism
            and self.donor_organism.corresponding_old_id
            != donor_organism.corresponding_old_id
        ):
            logging.log(
                logging.WARNING,
                f"WARNING: Specimen with ID {self.corresponding_old_id} already has "
                f"Donor Organism {self.donor_organism.corresponding_old_id} already "
                f"associated. It is being replaced.",
            )

        self.donor_organism = donor_organism

    def to_dictionary(self):
        dictionary_representation = deepcopy(self.__dict__)
        dictionary_representation[
            "donor_organism"
        ] = self.donor_organism.corresponding_old_id
        return dictionary_representation
