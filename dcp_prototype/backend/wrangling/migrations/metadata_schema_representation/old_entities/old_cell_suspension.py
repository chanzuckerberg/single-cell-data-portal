from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_specimen_from_organism import (
    OldSpecimenFromOrganism,
)

from dcp_prototype.backend.ledger.code.common.ledger_orm import BiosamplePrep
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.biosample_prep import (
    BiosamplePrep,
)
from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_transformer,
)
import logging


class OldCellSuspension:
    def __init__(self):
        self.corresponding_old_id = None

        self.specimen_from_organism = None

    def populate_from_dcp_one_json_data_frame(self, row):
        self.corresponding_old_id = row.get("provenance.document_id")

    def set_specimen_from_organism(
        self, specimen_from_organism: OldSpecimenFromOrganism
    ):
        if (
            self.specimen_from_organism
            and self.specimen_from_organism.corresponding_old_id
            != specimen_from_organism.corresponding_old_id
        ):
            logging.log(
                logging.WARNING,
                f"WARNING: Cell Suspension with ID {self.corresponding_old_id} already "
                f"has Specimen {self.specimen_from_organism.corresponding_old_id} "
                f"already associated. It is being replaced.",
            )

        self.specimen_from_organism = specimen_from_organism

    def to_dictionary(self):
        dictionary_representation = self.__dict__
        dictionary_representation[
            "specimen_from_organism"
        ] = self.specimen_from_organism.corresponding_old_id
        return dictionary_representation

    def convert_to_new_entity(self):
        biosample_prep_id = hca_accession_transformer(
            BiosamplePrep.__class__.__name__, self.corresponding_old_id
        )
        biosample_prep = BiosamplePrep(
            id=biosample_prep_id,
            category="primary_specimen",
            organ_ontology=self.specimen_from_organism.organ,
            developmental_stage=self.specimen_from_organism.donor_organism.developmental_stage,
            disease_ontology=self.specimen_from_organism.donor_organism.disease_onotology_label,
        )
        return biosample_prep
