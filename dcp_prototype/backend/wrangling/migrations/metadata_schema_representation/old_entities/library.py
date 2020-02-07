from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.library_prep_protocol import (
    LibraryPrepProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.sequencing_protocol import (
    SequencingProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.project import (
    Project,
)


class Library:
    def __init__(self):
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

    def is_complete(self):
        return self.library_prep_protocol and self.project and self.sequencing_protocol

    def to_dictionary(self):
        return {
            "library_prep_protocol": self.library_prep_protocol.corresponding_old_id,
            "project": self.project.corresponding_old_id,
            "sequencing_protocol": self.sequencing_protocol.corresponding_old_id,
        }

    def __eq__(self, other):
        return (
            self.library_prep_protocol.corresponding_old_id
            == other.library_prep_protocol.corresponding_old_id
            and self.project.corresponding_old_id == other.project.corresponding_old_id
            and self.sequencing_protocol.corresponding_old_id
            == other.sequencing_protocol.corresponding_old_id
        )
