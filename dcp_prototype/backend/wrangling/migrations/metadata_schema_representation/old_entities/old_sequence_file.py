from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_sequencing_protocol import (  # noqa
    OldSequencingProtocol,
)
from dcp_prototype.backend.wrangling.migrations.utils.id_generator import hca_accession_transformer
from dcp_prototype.backend.ledger.code.common.ledger_orm import SequenceFile
import logging
from copy import deepcopy


class OldSequenceFile:
    def __init__(self):
        self.corresponding_old_id = None
        self.filename = None
        self.file_format = None
        self.flowcell_id = None
        self.lane_index = None
        self.read_index = None
        self.s3_uri = None

        self.sequencing_protocol = None

    def populate_from_dcp_one_json_data_frame(self, row):
        self.corresponding_old_id = row.get("provenance.document_id")
        self.filename = row.get("file_core.file_name")
        self.file_format = row.get("file_core.format")
        self.lane_index = row.get("lane_index")
        self.read_index = row.get("read_index")

    def set_s3_uri(self, uri):
        self.s3_uri = uri + self.filename

    def set_sequencing_protocol(self, sequencing_protocol: OldSequencingProtocol):
        if (
            self.sequencing_protocol
            and self.sequencing_protocol.corresponding_old_id != sequencing_protocol.corresponding_old_id
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
        dictionary_representation = deepcopy(self.__dict__)
        dictionary_representation["sequencing_protocol"] = self.sequencing_protocol.corresponding_old_id
        return dictionary_representation

    def convert_to_new_entity(self):
        sequence_file_id = hca_accession_transformer(SequenceFile.__name__, self.corresponding_old_id)

        sequence_file = SequenceFile(
            id=sequence_file_id,
            filename=self.filename,
            file_format=self.file_format,
            flowcell_id=self.flowcell_id,
            lane_index=self.lane_index,
            read_index=self.read_index,
            s3_uri=self.s3_uri,
        )

        return sequence_file
