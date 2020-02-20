from dcp_prototype.backend.ledger.code.common.ledger_orm import AnalysisFile
from dcp_prototype.backend.wrangling.migrations.utils.id_generator import hca_accession_transformer
from copy import deepcopy


class OldAnalysisFile:
    def __init__(self):
        self.corresponding_old_id = None
        self.filename = None
        self.file_format = None
        self.s3_uri = None

    def populate_from_dcp_one_json_data_frame(self, row):
        self.corresponding_old_id = row.get("provenance.document_id")
        self.filename = row.get("file_core.file_name")
        self.file_format = row.get("file_core.format")

    def set_s3_uri(self, uri):
        self.s3_uri = uri + self.filename

    def to_dictionary(self):
        return deepcopy(self.__dict__)

    def convert_to_new_entity(self):
        analysis_file_id = hca_accession_transformer(AnalysisFile.__name__, self.corresponding_old_id)

        analysis_file = AnalysisFile(
            id=analysis_file_id, filename=self.filename, file_format=self.file_format, s3_uri=self.s3_uri,
        )

        return analysis_file
