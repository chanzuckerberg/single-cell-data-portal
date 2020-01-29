from utils.id_generator import hca_accession_generator


class Library:
    def __init__(self):
        self.id = hca_accession_generator(self.__class__.__name__)
