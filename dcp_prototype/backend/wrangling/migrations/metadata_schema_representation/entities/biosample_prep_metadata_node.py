from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_generator,
)
from dcp_prototype.backend.ledger.code.common.ledger_orm import BiosamplePrep

class BiosamplePrepMetadataNode(BiosamplePrep):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def to_dictionary(self):
        dictionary_rep = {}
        for key, value in self.__dict__.items():
            if key != "library_prep_protocols":
                dictionary_rep[key] = value
        return dictionary_rep
