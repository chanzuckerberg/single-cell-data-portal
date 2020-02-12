from dcp_prototype.backend.ledger.code.common.ledger_orm import SequencingProtocol


class SequencingProtocolMetadataNode(SequencingProtocol):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def to_dictionary(self):
        return self.__dict__
