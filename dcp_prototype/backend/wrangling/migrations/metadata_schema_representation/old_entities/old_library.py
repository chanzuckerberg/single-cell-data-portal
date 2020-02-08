from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_library_prep_protocol import (
    OldLibraryPrepProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_sequencing_protocol import (
    OldSequencingProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_project import (
    OldProject,
)
from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_transformer,
    hca_accession_generator,
)

from dcp_prototype.backend.ledger.code.common.ledger_orm import (
    Library,
    LibraryPrepProtocol,
    Project,
    SequencingProtocol,
)


class OldLibrary:
    def __init__(self):
        self.library_prep_protocol = None
        self.project = None
        self.sequencing_protocol = None

    def populate_associated_library_prep_protocol(
        self, library_prep_protocol: OldLibraryPrepProtocol
    ):
        self.library_prep_protocol = library_prep_protocol

    def populate_associated_project(self, project: OldProject):
        self.project = project

    def populate_associated_sequencing_protocol(
        self, sequencing_protocol: OldSequencingProtocol
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

    def convert_to_new_entity(self):
        library_id = hca_accession_generator(Library.__class__.__name__)
        library_prep_id = hca_accession_transformer(
            LibraryPrepProtocol.__class__.__name__,
            self.library_prep_protocol.corresponding_old_id,
        )
        project_id = hca_accession_transformer(
            Project.__class__.__name__, self.project.corresponding_old_id
        )
        sequencing_protocol_id = hca_accession_transformer(
            SequencingProtocol.__class__.__name__,
            self.sequencing_protocol.corresponding_old_id,
        )

        library = Library(
            id=library_id,
            library_prep_protocol_id=library_prep_id,
            sequencing_protocol_id=sequencing_protocol_id,
            project_id=project_id,
        )

        return library
