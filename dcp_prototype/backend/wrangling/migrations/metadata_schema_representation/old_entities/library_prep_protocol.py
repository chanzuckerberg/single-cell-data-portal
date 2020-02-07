from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.cell_suspension import (
    CellSuspension,
)
import logging


class LibraryPrepProtocol:
    def __init__(self):
        self.corresponding_old_id = None
        self.input_nucleic_acid_molecule_ontology = None
        self.library_construction_method_ontology = None
        self.nucleic_acid_source = None
        self.end_bias = None

        self.cell_suspension = None

    def populate_from_dcp_one_json_data_frame(self, row):
        self.input_nucleic_acid_molecule_ontology = row.get(
            "input_nucleic_acid_molecule_ontology.ontology"
        )
        self.library_construction_method_ontology = row.get(
            "library_construction_method.ontology"
        )
        self.nucleic_acid_source = row.get("nucleic_acid_source")
        self.end_bias = row.get("end_bias")
        self.corresponding_old_id = row.get("provenance.document_id")

    def set_cell_suspension(self, cell_suspension: CellSuspension):
        if (
            self.cell_suspension
            and self.cell_suspension.corresponding_old_id
            != cell_suspension.corresponding_old_id
        ):
            logging.log(
                logging.WARNING,
                f"WARNING: Library Prep with ID {self.corresponding_old_id} already "
                f"has Cell Suspension {self.cell_suspension.corresponding_old_id} "
                f"already associated. It is being replaced.",
            )

        self.cell_suspension = cell_suspension

    def to_dictionary(self):
        dictionary_representation = self.__dict__
        dictionary_representation[
            "cell_suspension"
        ] = self.cell_suspension.corresponding_old_id
        return dictionary_representation
