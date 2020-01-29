from utils.id_generator import hca_accession_generator


class BiosamplePrep:
    def __init__(
        self,
        category,
        organ_ontology,
        developmental_stage,
        disease_ontology_label,
    ):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.category = category
        self.organ_ontology = organ_ontology
        self.developmental_stage = developmental_stage
        self.disease_onotology_label = disease_ontology_label
