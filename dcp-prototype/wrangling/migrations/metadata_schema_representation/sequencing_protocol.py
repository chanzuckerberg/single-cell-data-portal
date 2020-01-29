from utils.id_generator import hca_accession_generator

class SequencingProtocol:
    def __init__(self, paird_end_sequencing, instrument_manufacturer_model):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.paired_end_sequencing = paird_end_sequencing
        self.instrument_manufacturer_model = instrument_manufacturer_model
