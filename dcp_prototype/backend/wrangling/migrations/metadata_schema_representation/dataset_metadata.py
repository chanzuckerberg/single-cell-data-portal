from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.biosample_prep import (
    BiosamplePrep,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.library_prep_protocol import (
    LibraryPrepProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.library import (
    Library,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.project import (
    Project,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.contributor import (
    Contributor,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.sequence_file import (
    SequenceFile,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.sequencing_protocol import (
    SequencingProtocol,
)

from dcp_prototype.backend.wrangling.migrations.utils.util import merge_dictionary_into
from pandas import DataFrame, ExcelWriter


class DatasetMetadata:
    def __init__(self):
        self.map_biosample_prep_ids_from_old_to_new = {}
        self.map_library_prep_protocol_ids_from_old_to_new = {}
        self.map_project_ids_from_old_to_new = {}
        self.map_contributor_ids_from_old_to_new = {}
        self.map_sequencing_protocol_ids_from_old_to_new = {}
        self.map_sequence_file_ids_from_old_to_new = {}
        self.libraries = []

    def export_to_spreadsheet(self):
        biosample_data_frame_dict = self.transform_metadata_entities_to_data_frame(
            self.map_biosample_prep_ids_from_old_to_new.values()
        )
        biosample_metadata_data_frame = DataFrame(
            biosample_data_frame_dict, columns=biosample_data_frame_dict.keys()
        )

        library_prep_data_frame_dict = self.transform_metadata_entities_to_data_frame(
            self.map_library_prep_protocol_ids_from_old_to_new.values()
        )
        library_prep_metadata_data_frame = DataFrame(
            library_prep_data_frame_dict, columns=library_prep_data_frame_dict.keys()
        )

        library_data_frame_dict = self.transform_metadata_entities_to_data_frame(
            self.libraries
        )
        library_metadata_data_frame = DataFrame(
            library_data_frame_dict, columns=library_data_frame_dict.keys()
        )

        project_data_frame_dict = self.transform_metadata_entities_to_data_frame(
            self.map_project_ids_from_old_to_new.values()
        )
        project_metadata_data_frame = DataFrame(
            project_data_frame_dict, columns=project_data_frame_dict.keys()
        )

        contributor_data_frame_dict = self.transform_metadata_entities_to_data_frame(
            self.map_contributor_ids_from_old_to_new.values()
        )
        contributor_metadata_data_frame = DataFrame(
            contributor_data_frame_dict, columns=contributor_data_frame_dict.keys()
        )

        sequencing_protocol_data_frame_dict = self.transform_metadata_entities_to_data_frame(
            self.map_sequencing_protocol_ids_from_old_to_new.values()
        )
        sequencing_protocol_metadata_data_frame = DataFrame(
            sequencing_protocol_data_frame_dict,
            columns=sequencing_protocol_data_frame_dict.keys(),
        )

        sequence_file_data_frame_dict = self.transform_metadata_entities_to_data_frame(
            self.map_sequence_file_ids_from_old_to_new.values()
        )
        sequence_file_metadata_data_frame = DataFrame(
            sequence_file_data_frame_dict, columns=sequence_file_data_frame_dict.keys()
        )

        with ExcelWriter("metadata_spreadsheet.xlsx") as excel_writer:
            biosample_metadata_data_frame.to_excel(
                excel_writer, sheet_name="Biosample Prep", index=False
            )
            library_prep_metadata_data_frame.to_excel(
                excel_writer, sheet_name="Library Prep Protocol", index=False
            )
            library_metadata_data_frame.to_excel(
                excel_writer, sheet_name="Library", index=False
            )
            project_metadata_data_frame.to_excel(
                excel_writer, sheet_name="Project", index=False
            )
            contributor_metadata_data_frame.to_excel(
                excel_writer, sheet_name="Contributor", index=False
            )
            sequencing_protocol_metadata_data_frame.to_excel(
                excel_writer, sheet_name="Sequencing Protocol", index=False
            )
            sequence_file_metadata_data_frame.to_excel(
                excel_writer, sheet_name="Sequence File", index=False
            )

    def transform_metadata_entities_to_data_frame(self, metadata_entities_list):
        dict_form_of_data = {}
        for metadata_entity in metadata_entities_list:
            dict_form_of_data = merge_dictionary_into(
                dict_form_of_data, metadata_entity.to_dictionary()
            )
        return dict_form_of_data

    def parse_row_of_metadata(self, row):
        biosample_prep_old_id = row.get(
            "donor_organism.biomaterial_core.biomaterial_name"
        )
        if (
            biosample_prep_old_id
            not in self.map_biosample_prep_ids_from_old_to_new.keys()
        ):
            biosample_prep = BiosamplePrep()
            biosample_prep.populate_from_dcp_one_data_row(row)
            self.map_biosample_prep_ids_from_old_to_new[
                biosample_prep_old_id
            ] = biosample_prep

        library_prep_protocol_old_id = row.get(
            "library_preparation_protocol.protocol_core.protocol_name"
        )
        if (
            library_prep_protocol_old_id
            not in self.map_library_prep_protocol_ids_from_old_to_new.keys()
        ):
            library_prep_protocol = LibraryPrepProtocol()
            library_prep_protocol.populate_from_dcp_one_data_row(row)

            library_prep_protocol.populate_associated_biosample(
                self.map_biosample_prep_ids_from_old_to_new[biosample_prep_old_id]
            )

            self.map_library_prep_protocol_ids_from_old_to_new[
                library_prep_protocol_old_id
            ] = library_prep_protocol

        contact_old_id = row.get("project.contributors.contact_name")
        if contact_old_id not in self.map_contributor_ids_from_old_to_new.keys():
            contributor = Contributor()
            contributor.populate_from_dcp_one_data_row(row)
            self.map_contributor_ids_from_old_to_new[contact_old_id] = contributor

        project_old_id = row.get("project.project_core.project_title")
        if project_old_id not in self.map_project_ids_from_old_to_new.keys():
            project = Project()
            project.populate_from_dcp_one_data_row(row)

            project.add_associated_contributor(
                self.map_contributor_ids_from_old_to_new[contact_old_id]
            )

            self.map_project_ids_from_old_to_new[project_old_id] = project

        sequencing_file_old_id = row.get("*.file_core.file_name")
        if (
            sequencing_file_old_id
            not in self.map_sequence_file_ids_from_old_to_new.keys()
        ):
            sequencing_file = SequenceFile()
            sequencing_file.populate_from_dcp_one_data_row(row)

            self.map_sequence_file_ids_from_old_to_new[
                sequencing_file_old_id
            ] = sequencing_file

        sequencing_protocol_old_id = row.get(
            "sequencing_protocol.protocol_core.protocol_name"
        )
        if (
            sequencing_protocol_old_id
            not in self.map_sequencing_protocol_ids_from_old_to_new.keys()
        ):
            sequencing_protocol = SequencingProtocol()
            sequencing_protocol.populate_from_dcp_one_data_row(row)

            sequencing_protocol.add_associated_sequencing_file(
                self.map_sequence_file_ids_from_old_to_new[sequencing_file_old_id]
            )

            self.map_sequencing_protocol_ids_from_old_to_new[
                sequencing_protocol_old_id
            ] = sequencing_protocol

        library = Library()
        library.populate_associated_project(
            self.map_project_ids_from_old_to_new[project_old_id]
        )
        library.populate_associated_sequencing_protocol(
            self.map_sequencing_protocol_ids_from_old_to_new[sequencing_protocol_old_id]
        )
        library.populate_associated_library_prep_protocol(
            self.map_library_prep_protocol_ids_from_old_to_new[
                library_prep_protocol_old_id
            ]
        )
        if library not in self.libraries:
            self.libraries.append(library)
