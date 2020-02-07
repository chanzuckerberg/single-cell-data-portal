from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.sequencing_protocol import (
    SequencingProtocol,
)
import logging


class SequenceFile:
    def __init__(self):
        self.corresponding_old_id = None
        self.filename = None
        self.file_format = None
        self.flowcell_id = None
        self.lane_index = None
        self.read_index = None

        self.sequencing_protocol = None

    def populate_from_dcp_one_json_data_frame(self, row):
        self.corresponding_old_id = row.get("provenance.document_id")
        self.filename = row.get("file_core.file_name")
        self.file_format = row.get("file_core.format")
        self.lane_index = row.get("lane_index")
        self.read_index = row.get("read_index")

    def set_sequencing_protocol(self, sequencing_protocol: SequencingProtocol):
        if (
            self.sequencing_protocol
            and self.sequencing_protocol.corresponding_old_id
            != sequencing_protocol.corresponding_old_id
        ):
            logging.log(
                logging.WARNING,
                f"WARNING: Sequence File with ID {self.corresponding_old_id} already "
                f"has Sequencing Protocol "
                f"{self.sequencing_protocol.corresponding_old_id} already associated. "
                f"It is being replaced.",
            )

        self.sequencing_protocol = sequencing_protocol

    def to_dictionary(self):
        dictionary_representation = self.__dict__
        dictionary_representation[
            "sequencing_protocol"
        ] = self.sequencing_protocol.corresponding_old_id
        return dictionary_representation
