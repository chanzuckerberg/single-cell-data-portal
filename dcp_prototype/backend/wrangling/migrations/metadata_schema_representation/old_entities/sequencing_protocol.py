class SequencingProtocol:
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
        return self.__dict__