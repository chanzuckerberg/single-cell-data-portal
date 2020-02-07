from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_generator,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.biosample_prep import (  # noqa
    BiosamplePrep,
)
import logging


class LibraryPrepProtocol:
    def __init__(self):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.input_nucleic_acid_molecule_ontology = None
        self.library_construction_method_ontology = None
        self.nucleic_acid_source = None
        self.end_bias = None

        self.biosample_prep = None

    def populate_from_dcp_one_data_row(self, row):
        self.input_nucleic_acid_molecule_ontology = row.get(
            "library_preparation_protocol.input_nucleic_acid_molecule_ontology."
            "ontology_label"
        )
        self.library_construction_method_ontology = row.get(
            "library_preparation_protocol.library_construction_method.ontology_label"
        )
        self.nucleic_acid_source = row.get(
            "library_preparation_protocol.nucleic_acid_source"
        )
        self.end_bias = row.get("library_preparation_protocol.end_bias")

    def populate_from_dcp_one_json_data_frame(self, row):
        self.input_nucleic_acid_molecule_ontology = row.get(
            "input_nucleic_acid_molecule_ontology.ontology"
        )
        self.library_construction_method_ontology = row.get(
            "library_construction_method.ontology"
        )
        self.nucleic_acid_source = row.get("nucleic_acid_source")
        self.end_bias = row.get("end_bias")

    def populate_associated_biosample(self, biosample_prep: BiosamplePrep):
        if self.biosample_prep and self.biosample_prep.id != biosample_prep.id:
            # Log warning if there are multiple biosample preps associated with this
            # library prep protocol
            logging.log(
                logging.WARNING,
                f"WARNING: Library prep {self.id} is already associated with biosample "
                f"prep {self.biosample_prep.id}. There are multiple biosample preps "
                f"associated and now it is being replaced by biosample prep with id "
                f"{biosample_prep.id}!",
            )
        self.biosample_prep = biosample_prep

    def to_dictionary(self):
        dictionary_representation = self.__dict__
        dictionary_representation["biosample_prep"] = self.biosample_prep.id
        return dictionary_representation
