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

    def populate_from_dcp_one_data_row(self, row):
        self.organ_ontology = row.get("specimen_from_organism.organ.ontology_label")
        self.developmental_stage = row.get(
            "donor_organism.development_stage.ontology_label"
        )
        self.disease_onotology_label = row.get("donor_organism.diseases.ontology_label")

    def to_dictionary(self):
        return self.__dict__
