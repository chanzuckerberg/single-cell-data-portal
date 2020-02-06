from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_generator,
)


class BiosamplePrep:
    def __init__(self):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.category = "primary_specimen"
        self.organ_ontology = None
        self.developmental_stage = None
        self.disease_onotology_label = None
        self.other_associated_ids = []

    def add_associated_id(self, id):
        self.other_associated_ids.append(id)

    def populate_from_dcp_one_data_row(self, row):
        self.organ_ontology = row.get("specimen_from_organism.organ.ontology_label")
        self.developmental_stage = row.get(
            "donor_organism.development_stage.ontology_label"
        )
        self.disease_onotology_label = row.get("donor_organism.diseases.ontology_label")

    def populate_from_dcp_one_json_data_frame(self, row):
        self.organ_ontology = None  # TODO
        self.developmental_stage = row.get("development_stage.ontology")
        self.disease_onotology_label = row.get("diseases.0.ontology")

    def to_dictionary(self):
        dictionary_rep = {}
        for key, value in self.__dict__.items():
            if key != "other_associated_ids":
                dictionary_rep[key] = value
        return dictionary_rep
