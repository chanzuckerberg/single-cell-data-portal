from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.specimen_from_organism import (
    SpecimenFromOrganism,
)
import logging


class CellSuspension:
    def __init__(self):
        self.corresponding_old_id = None

        self.specimen_from_organism = None

    def populate_from_dcp_one_json_data_frame(self, row):
        self.corresponding_old_id = row.get("provenance.document_id")

    def set_specimen_from_organism(self, specimen_from_organism: SpecimenFromOrganism):
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
