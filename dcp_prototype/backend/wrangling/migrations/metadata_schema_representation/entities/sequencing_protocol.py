from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_generator,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.sequence_file import (  # noqa
    SequenceFile,
)
from dcp_prototype.backend.wrangling.migrations.utils.util import merge_dictionary_into


class SequencingProtocol:
    def __init__(self):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.paired_end_sequencing = None
        self.instrument_manufacturer_model = None

        self.sequence_files = []

    def populate_from_dcp_one_data_row(self, row):
        self.paired_end_sequencing = row.get("sequencing_protocol.paired_end")
        self.instrument_manufacturer_model = row.get(
            "sequencing_protocol.instrument_manufacturer_model.ontology_label"
        )

    def __eq__(self, other):
        return (
            self.paired_end_sequencing == other.paired_end_sequencing
            and self.instrument_manufacturer_model
            == other.instrument_manufacturer_model
        )

    def add_associated_sequencing_file(self, sequence_file: SequenceFile):
        self.sequence_files.append(sequence_file)

    def to_dictionary(self):
        dictionary_representation = {}
        for sequence_file in self.sequence_files:
            for key, value in self.__dict__.items():
                if key == "sequence_files":
                    merge_dictionary_into(
                        dictionary_representation, {key: sequence_file.id}
                    )
                else:
                    merge_dictionary_into(dictionary_representation, {key: value})
        return dictionary_representation
