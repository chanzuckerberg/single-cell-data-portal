class OldDonorOrganism:
    def __init__(self):
        self.corresponding_old_id = None
        self.developmental_stage = None
        self.disease_onotology_label = None

    def populate_from_dcp_one_json_data_frame(self, row):
        self.corresponding_old_id = row.get("provenance.document_id")
        self.developmental_stage = row.get("development_stage.ontology")
        self.disease_onotology_label = row.get("diseases.0.ontology")

    def to_dictionary(self):
        return self.__dict__
