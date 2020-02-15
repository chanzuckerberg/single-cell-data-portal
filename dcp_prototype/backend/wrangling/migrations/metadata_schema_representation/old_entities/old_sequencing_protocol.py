from dcp_prototype.backend.ledger.code.common.ledger_orm import SequencingProtocol
from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_transformer,
)
from copy import deepcopy


class OldSequencingProtocol:
    def __init__(self):
        self.corresponding_old_id = None
        self.paired_end_sequencing = None
        self.instrument_manufacturer_model = None

    def populate_from_dcp_one_json_data_frame(self, row):
        self.corresponding_old_id = row.get("provenance.document_id")
        self.paired_end_sequencing = row.get("paired_end")
        self.instrument_manufacturer_model = row.get(
            "instrument_manufacturer_model.ontology"
        )

    def to_dictionary(self):
        return deepcopy(self.__dict__)

    def convert_to_new_entity(self):
        protocol_id = hca_accession_transformer(
            SequencingProtocol.__name__, self.corresponding_old_id
        )

        sequencing_protocol = SequencingProtocol(
            id=protocol_id,
            paired_end_sequencing=self.paired_end_sequencing,
            instrument_manufacturer_model=self.instrument_manufacturer_model,
        )

        return sequencing_protocol
