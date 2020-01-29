from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_generator,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.project import (
    Project,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.library_prep_protocol import (
    LibraryPrepProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.sequencing_protocol import (
    SequencingProtocol,
)


class Library:
    def __init__(self):
        self.id = hca_accession_generator(self.__class__.__name__)
        self.library_prep_protocol = None
        self.project = None
        self.sequencing_protocol = None

    def populate_associated_library_prep_protocol(
        self, library_prep_protocol: LibraryPrepProtocol
    ):
        self.library_prep_protocol = library_prep_protocol

    def populate_associated_project(self, project: Project):
        self.project = project

    def populate_associated_sequencing_protocol(
        self, sequencing_protocol: SequencingProtocol
    ):
        self.sequencing_protocol = sequencing_protocol

    def __eq__(self, other):
        return (
            self.library_prep_protocol == other.library_prep_protocol
            and self.project == other.project
            and self.sequencing_protocol == other.sequencing_protocol
        )

    def pretty_print(self):
        print(self.library_prep_protocol.pretty_print())
        print(self.project.pretty_print())
        print(self.sequencing_protocol.pretty_print())

    def to_dictionary(self):
        return {
            "id": self.id,
            "library_prep_protocol": self.library_prep_protocol.id,
            "project": self.project.id,
            "sequencing_protocol": self.sequencing_protocol.id,
        }
