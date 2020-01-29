from utils.id_generator import hca_accession_generator


class SequenceFile:
    def __init__(self, filename, file_format, flowcell_id, lane_index, read_index):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.filename = filename
        self.file_format = file_format
        self.flowcell_id = flowcell_id
        self.lane_index = lane_index
        self.read_index = read_index
