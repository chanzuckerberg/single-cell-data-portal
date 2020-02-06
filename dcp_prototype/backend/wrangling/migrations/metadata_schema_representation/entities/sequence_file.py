from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_generator,
)


class SequenceFile:
    def __init__(self):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.filename = None
        self.file_format = None
        self.flowcell_id = None
        self.lane_index = None
        self.read_index = None
        self.s3_uri = None

    def populate_from_dcp_one_data_row(self, row):
        self.filename = row.get("*.file_core.file_name")
        self.file_format = row.get("sequence_file.file_core.format")
        self.lane_index = row.get("sequence_file.lane_index")
        self.read_index = row.get("sequence_file.read_index")

    def set_s3_uri(self, s3_uri):
        self.s3_uri = s3_uri

    def populate_from_dcp_one_json_data_frame(self, row):
        self.filename = row.get("file_core.file_name")
        self.file_format = row.get("file_core.format")
        self.lane_index = row.get("lane_index")
        self.read_index = row.get("read_index")

    def to_dictionary(self):
        return self.__dict__
