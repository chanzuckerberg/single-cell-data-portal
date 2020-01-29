from utils.id_generator import hca_accession_generator


class LibraryPrepProtocol:
    def __init__(
        self,
        input_nucleic_acid_molecule,
        library_construction_method_ontology,
        nucleic_acid_source,
        end_bias,
    ):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.input_nucleic_acid_molecule = input_nucleic_acid_molecule
        self.library_construction_method_ontology = library_construction_method_ontology
        self.nucleic_acid_source = nucleic_acid_source
        self.end_bias = end_bias
