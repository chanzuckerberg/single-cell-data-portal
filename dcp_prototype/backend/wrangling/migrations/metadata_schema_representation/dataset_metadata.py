from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.biosample_prep import (  # noqa
    BiosamplePrep,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.library_prep_protocol import (  # noqa
    LibraryPrepProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.library import (  # noqa
    Library,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.project import (  # noqa
    Project,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.contributor import (  # noqa
    Contributor,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.sequence_file import (  # noqa
    SequenceFile,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.entities.sequencing_protocol import (  # noqa
    SequencingProtocol,
)

from dcp_prototype.backend.wrangling.migrations.utils.util import merge_dictionary_into
from pandas import DataFrame, ExcelWriter
import random


class DatasetMetadata:
    """
    DatasetMetadata captures the entire experimental metadata of a single project and
    links the entities that make up the experimental metadata graph together. It
    provides a capability to also export the metadata structure into a spreadsheet.

    The object has two modes: editing mode and publishing mode. In editing mode, the
    structure is still being built up and the linking between the entities is being
    established on the fly. In publish mode, the structure is able to be exported to a
    spreadsheet. We separate the two modes because publishing metadata that has not yet
    been linked will lead to confusion about the accuracy of the spreadsheet.
    """

    def __init__(self):
        self.map_biosample_prep_ids_from_old_to_new = {}
        self.map_library_prep_protocol_ids_from_old_to_new = {}
        self.map_project_ids_from_old_to_new = {}
        self.map_contributor_ids_from_old_to_new = {}
        self.map_sequencing_protocol_ids_from_old_to_new = {}
        self.map_sequence_file_ids_from_old_to_new = {}
        self.libraries = []

        self.saved_full_metadata_structure = {}

        self.publish_mode = False

    def export_to_spreadsheet(self, spreadsheet_filename):
        """
        Write the experimental metadata representation graph to the given spreadsheet
        filename where each sheet represents a single type of entity and each row is a
        unique instance of that entity that exists in the experimental graph.
        """

        if not self.publish_mode:
            raise Exception(
                "Please run the save() method on the DataMetadata instance so that the "
                "linking can be finalized prior to export to spreadsheet."
            )
        with ExcelWriter(spreadsheet_filename) as excel_writer:
            for (
                sheet_name,
                entity_entries,
            ) in self.saved_full_metadata_structure.items():
                entity_transformation_to_dictionary_format = self._transform_metadata_entities_to_data_frame(  # noqa
                    entity_entries
                )
                entity_transformation_to_data_frame_format = DataFrame(
                    entity_transformation_to_dictionary_format,
                    columns=entity_transformation_to_dictionary_format.keys(),
                )
                entity_transformation_to_data_frame_format.to_excel(
                    excel_writer, sheet_name=sheet_name, index=False
                )

    def save(self):
        """
        Extracts only the entities themselves from the maps that have been keeping track
        of the mapping between the inputted metadata entities and the ones that are now
        being created as abstract classes. Since the spreadsheet only requires the
        entities and the entities self-sufficiently capture the edges/relationships
        between the entities, flip the state of the object to publish mode so that it is
        clear that the entities being outputted reflect the full contents of a metadata
        spreadsheet.
        """

        self.saved_full_metadata_structure[
            "Biosample Prep"
        ] = self.map_biosample_prep_ids_from_old_to_new.values()
        self.saved_full_metadata_structure[
            "Library Prep Protocol"
        ] = self.map_library_prep_protocol_ids_from_old_to_new.values()
        self.saved_full_metadata_structure["Library"] = self.libraries
        self.saved_full_metadata_structure[
            "Project"
        ] = self.map_project_ids_from_old_to_new.values()
        self.saved_full_metadata_structure[
            "Contributors"
        ] = self.map_contributor_ids_from_old_to_new.values()
        self.saved_full_metadata_structure[
            "Sequencing Protocol"
        ] = self.map_sequencing_protocol_ids_from_old_to_new.values()
        self.saved_full_metadata_structure[
            "Sequence File"
        ] = self.map_sequence_file_ids_from_old_to_new.values()

        if self.saved_full_metadata_structure["Project"]:
            self.project_name = list(self.saved_full_metadata_structure["Project"])[
                0
            ].project_short_name
        else:
            self.project_name = "For some reason this project got no name."

        self.publish_mode = True

    def parse_flattened_row_of_json(self, row, entity_type):
        if entity_type is "biosample_prep":
            biosample_prep_old_id = row.get("provenance.document_id")
            if (
                biosample_prep_old_id
                not in self.map_biosample_prep_ids_from_old_to_new.keys()
            ):
                biosample_prep = BiosamplePrep()
                biosample_prep.populate_from_dcp_one_json_data_frame(row)
                self.map_biosample_prep_ids_from_old_to_new[
                    biosample_prep_old_id
                ] = biosample_prep

        if entity_type is "specimen":
            for key, value in self.map_biosample_prep_ids_from_old_to_new.items():
                value.organ_ontology = row.get("organ.ontology")
                self.map_biosample_prep_ids_from_old_to_new[key] = value

        if entity_type is "other_biosample":
            for key, value in self.map_biosample_prep_ids_from_old_to_new.items():
                value.add_associated_id(row.get("provenance.document_id"))
                self.map_biosample_prep_ids_from_old_to_new[key] = value

        if entity_type is "library_preparation_protocol":
            library_prep_protocol_old_id = row.get("provenance.document_id")
            if (
                library_prep_protocol_old_id
                not in self.map_library_prep_protocol_ids_from_old_to_new.keys()
            ):
                library_prep_protocol = LibraryPrepProtocol()
                library_prep_protocol.populate_from_dcp_one_json_data_frame(row)
                self.map_library_prep_protocol_ids_from_old_to_new[
                    library_prep_protocol_old_id
                ] = library_prep_protocol

        if entity_type is "project":
            project_old_id = row.get("provenance.document_id")
            if project_old_id not in self.map_project_ids_from_old_to_new.keys():
                project = Project()
                project.populate_from_dcp_one_json_data_frame(row)

                self.map_project_ids_from_old_to_new[project_old_id] = project

                self.map_project_ids_from_old_to_new[project_old_id] = project

            for contributor in project.contributors:
                self.map_contributor_ids_from_old_to_new[contributor.id] = contributor

        if entity_type is "sequence_file":
            sequencing_file_old_id = row.get("provenance.document_id")
            if (
                sequencing_file_old_id
                not in self.map_sequence_file_ids_from_old_to_new.keys()
            ):
                sequencing_file = SequenceFile()
                sequencing_file.populate_from_dcp_one_json_data_frame(row)

                self.map_sequence_file_ids_from_old_to_new[
                    sequencing_file_old_id
                ] = sequencing_file

        if entity_type is "sequencing_protocol":
            print("IM HERE!")
            sequencing_protocol_old_id = row.get("provenance.document_id")
            if (
                sequencing_protocol_old_id
                not in self.map_sequencing_protocol_ids_from_old_to_new.keys()
            ):
                sequencing_protocol = SequencingProtocol()
                sequencing_protocol.populate_from_dcp_one_json_data_frame(row)

                self.map_sequencing_protocol_ids_from_old_to_new[
                    sequencing_protocol_old_id
                ] = sequencing_protocol
                print(self.map_sequencing_protocol_ids_from_old_to_new)

        if entity_type is "links":
            link_index = 0
            while row.get(f"links.{str(link_index)}.process"):
                link_prefix = f"links.{str(link_index)}."

                library = Library()
                library.populate_associated_project(
                    random.choice(list(self.map_project_ids_from_old_to_new.values()))
                )
                if row.get(f"{link_prefix}input_type") == "biomaterial":
                    protocol_index = 0
                    while row.get(
                        f"{link_prefix}protocols.{str(protocol_index)}.protocol_type"
                    ):
                        protocol_prefix = (
                            f"{link_prefix}protocols.{str(protocol_index)}."
                        )
                        protocol_type_key = f"{protocol_prefix}protocol_type"
                        protocol_id_key = f"{protocol_prefix}protocol_id"

                        if row.get(protocol_type_key) == "library_preparation_protocol":
                            library_prep_protocol_old_id = row.get(protocol_id_key)

                            library.populate_associated_library_prep_protocol(
                                self.map_library_prep_protocol_ids_from_old_to_new[
                                    library_prep_protocol_old_id
                                ]
                            )

                            biomaterial_index = 0
                            while row.get(
                                f"{link_prefix}inputs.{str(biomaterial_index)}"
                            ):
                                biosample_prep_old_id = row.get(
                                    f"{link_prefix}inputs.{str(biomaterial_index)}"
                                )

                                library_prep = self.map_library_prep_protocol_ids_from_old_to_new[
                                    library_prep_protocol_old_id
                                ]

                                for (
                                    old_id,
                                    biosample,
                                ) in (
                                    self.map_biosample_prep_ids_from_old_to_new.items()
                                ):
                                    if (
                                        biosample_prep_old_id
                                        in biosample.other_associated_ids
                                    ):
                                        library_prep.populate_associated_biosample(
                                            biosample
                                        )

                                biomaterial_index += 1

                        if row.get(protocol_type_key) == "sequencing_protocol":
                            sequencing_protocol_old_id = row.get(protocol_id_key)

                            library.populate_associated_sequencing_protocol(
                                self.map_sequencing_protocol_ids_from_old_to_new[
                                    sequencing_protocol_old_id
                                ]
                            )

                            if row.get(f"{link_prefix}output_type") == "file":
                                file_index = 0
                                while row.get(
                                    f"{link_prefix}outputs.{str(file_index)}"
                                ):
                                    sequencing_file_old_id = row.get(
                                        f"{link_prefix}outputs.{str(file_index)}"
                                    )

                                    sequencing_protocol = self.map_sequencing_protocol_ids_from_old_to_new[
                                        sequencing_protocol_old_id
                                    ]
                                    sequencing_protocol.add_associated_sequencing_file(
                                        self.map_sequence_file_ids_from_old_to_new[
                                            sequencing_file_old_id
                                        ]
                                    )
                                    print(self.map_sequencing_protocol_ids_from_old_to_new)

                                    file_index += 1

                        protocol_index += 1

                if library.is_complete() and library not in self.libraries:
                    self.libraries.append(library)

                link_index += 1

    def parse_row_of_metadata(self, row):
        """
        Read in a row of data from a metadata TSV file. Each row should represent the
        entire metadata for a single file in the dataset. This knowledge is used to
        generate the entities in the experimental graph that generated that particular
        file of data and stores it (as well as the relationships between the entities).
        """

        self.publish_mode = False

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

    def _transform_metadata_entities_to_data_frame(self, metadata_entities_list):
        """
        Transform each entity into its dictionary format and then merge them into a
        single dictionary where each key is an entity's attribute and the values are all
        the values of that attributes from all the entities in the given list.
        """

        dict_form_of_data = {}
        for metadata_entity in metadata_entities_list:
            dict_form_of_data = merge_dictionary_into(
                dict_form_of_data, metadata_entity.to_dictionary()
            )
        return dict_form_of_data
