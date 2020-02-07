from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.project import (
    Project,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.sequencing_protocol import (
    SequencingProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.library_prep_protocol import (
    LibraryPrepProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.library import (
    Library,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.cell_suspension import (
    CellSuspension,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.contributor import (
    Contributor,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.sequence_file import (
    SequenceFile,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.donor_organism import (
    DonorOrganism,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.specimen_from_organism import (
    SpecimenFromOrganism,
)


from dcp_prototype.backend.wrangling.migrations.utils.util import merge_dictionary_into
from pandas import DataFrame, ExcelWriter
import random


class DatasetMetadata:
    def __init__(self):
        self.donor_organisms = {}
        self.specimens = {}
        self.cell_suspensions = {}
        self.library_preps = {}
        self.libraries = []
        self.projects = {}
        self.contributors = {}
        self.files = {}
        self.sequencing_protocols = {}

        self.publish_mode = False
        self.old_metadata_structure = {}
        self.new_metadata_structure = {}

        self.project_name = ""

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
            for (sheet_name, entity_entries,) in self.old_metadata_structure.items():
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

        self.publish_mode = True

        self.old_metadata_structure["Donor Organism"] = self.donor_organisms.values()
        self.old_metadata_structure["Specimen From Organism"] = self.specimens.values()
        self.old_metadata_structure["Cell Suspensions"] = self.cell_suspensions.values()
        self.old_metadata_structure[
            "Library Prep Protocols"
        ] = self.library_preps.values()
        self.old_metadata_structure["Libraries"] = self.libraries
        self.old_metadata_structure["Projects"] = self.projects.values()
        self.old_metadata_structure["Contributors"] = self.contributors.values()
        self.old_metadata_structure[
            "Sequencing Protocols"
        ] = self.sequencing_protocols.values()
        self.old_metadata_structure["Sequence Files"] = self.files.values()

        if self.old_metadata_structure["Projects"]:
            self.project_name = list(self.old_metadata_structure["Projects"])[
                0
            ].project_short_name
        else:
            self.project_name = "For some reason this project got no name."

        self.publish_mode = True

    def parse_flattened_row_of_json(self, row, entity_type):
        if entity_type == "donor_organism":
            donor_organism = DonorOrganism()
            donor_organism.populate_from_dcp_one_json_data_frame(row)
            if donor_organism.corresponding_old_id not in self.donor_organisms.keys():
                self.donor_organisms[
                    donor_organism.corresponding_old_id
                ] = donor_organism

        if entity_type == "specimen_from_organism":
            specimen_from_organism = SpecimenFromOrganism()
            specimen_from_organism.populate_from_dcp_one_json_data_frame(row)
            if specimen_from_organism.corresponding_old_id not in self.specimens.keys():
                self.specimens[
                    specimen_from_organism.corresponding_old_id
                ] = specimen_from_organism

        if entity_type == "cell_suspension":
            cell_suspension = CellSuspension()
            cell_suspension.populate_from_dcp_one_json_data_frame(row)
            if cell_suspension.corresponding_old_id not in self.cell_suspensions.keys():
                self.cell_suspensions[
                    cell_suspension.corresponding_old_id
                ] = cell_suspension

        if entity_type == "library_preparation_protocol":
            library_prep = LibraryPrepProtocol()
            library_prep.populate_from_dcp_one_json_data_frame(row)
            if library_prep.corresponding_old_id not in self.library_preps.keys():
                self.library_preps[library_prep.corresponding_old_id] = library_prep

        if entity_type == "project":
            project = Project()
            project.populate_from_dcp_one_json_data_frame(row)
            if project.corresponding_old_id not in self.projects.keys():
                self.projects[project.corresponding_old_id] = project
            for contributor in project.contributors:
                if contributor.name not in self.contributors.keys():
                    self.contributors[contributor.name] = contributor

        if entity_type == "sequencing_protocol":
            protocol = SequencingProtocol()
            protocol.populate_from_dcp_one_json_data_frame(row)
            if protocol.corresponding_old_id not in self.sequencing_protocols.keys():
                self.sequencing_protocols[protocol.corresponding_old_id] = protocol

        if entity_type == "sequence_file":
            sequence_file = SequenceFile()
            sequence_file.populate_from_dcp_one_json_data_frame(row)
            if sequence_file.corresponding_old_id not in self.files.keys():
                self.files[sequence_file.corresponding_old_id] = sequence_file

        if entity_type == "links":
            link_index = 0

            while row.get(f"links.{str(link_index)}.process"):
                link_index_prefix = f"links.{str(link_index)}."

                if row.get(f"{link_index_prefix}input_type") == "file":
                    link_index += 1
                    continue

                attempted_new_library = Library()
                attempted_new_library.populate_associated_project(
                    random.choice(list(self.projects.values()))
                )

                if (
                    row.get(f"{link_index_prefix}input_type") == "biomaterial"
                    and row.get(f"{link_index_prefix}output_type") == "file"
                ):
                    protocol_index = 0
                    while row.get(
                        f"{link_index_prefix}protocols.{str(protocol_index)}.protocol_type"
                    ):
                        protocol_type = row.get(
                            f"{link_index_prefix}protocols.{str(protocol_index)}.protocol_type"
                        )
                        if protocol_type == "sequencing_protocol":
                            protocol_id = row.get(
                                f"{link_index_prefix}protocols.{str(protocol_index)}.protocol_id"
                            )
                            associated_sequence_protocol = self.sequencing_protocols.get(
                                protocol_id
                            )

                            # Add to library generation
                            attempted_new_library.populate_associated_sequencing_protocol(
                                associated_sequence_protocol
                            )

                            # Link together sequence files and sequence protocol
                            file_index = 0
                            while row.get(
                                f"{link_index_prefix}outputs.{str(file_index)}"
                            ):
                                sequence_file_id = row.get(
                                    f"{link_index_prefix}outputs.{str(file_index)}"
                                )

                                file = self.files.get(sequence_file_id)
                                if file:
                                    file.set_sequencing_protocol(
                                        associated_sequence_protocol
                                    )
                                    self.files[file.corresponding_old_id] = file

                                file_index += 1

                        if protocol_type == "library_preparation_protocol":
                            protocol_id = row.get(
                                f"{link_index_prefix}protocols.{str(protocol_index)}.protocol_id"
                            )

                            associated_library_prep = self.library_preps.get(
                                protocol_id
                            )

                            # Add to library generation
                            attempted_new_library.populate_associated_library_prep_protocol(
                                associated_library_prep
                            )

                            # Link cell suspension and library prep together
                            input_index = 0
                            while row.get(
                                f"{link_index_prefix}inputs.{str(input_index)}"
                            ):
                                suspension_id = row.get(
                                    f"{link_index_prefix}inputs.{str(input_index)}"
                                )
                                associated_library_prep.set_cell_suspension(
                                    self.cell_suspensions.get(suspension_id)
                                )

                                input_index += 1

                        protocol_index += 1

                if (
                    row.get(f"{link_index_prefix}input_type") == "biomaterial"
                    and row.get(f"{link_index_prefix}output_type") == "biomaterial"
                ):
                    # If a link has no protocols associated, then the link is between
                    # donor organism and specimen
                    if row.get(
                        f"{link_index_prefix}protocols.0.protocol_type"
                    ) == None or not row.get(
                        f"{link_index_prefix}protocols.0.protocol_type"
                    ):
                        output_index = 0
                        while row.get(
                            f"{link_index_prefix}outputs.{str(output_index)}"
                        ):
                            specimen_id = row.get(
                                f"{link_index_prefix}outputs.{str(output_index)}"
                            )
                            specimen = self.specimens.get(specimen_id)

                            # Link together specimen and donor organism
                            input_index = 0
                            while row.get(
                                f"{link_index_prefix}inputs.{str(input_index)}"
                            ):
                                donor_id = row.get(
                                    f"{link_index_prefix}inputs.{str(input_index)}"
                                )
                                specimen.set_donor_organism(
                                    self.donor_organisms.get(donor_id)
                                )

                                input_index += 1

                            output_index += 1

                    # If a link has protocols associated with it, then the link is
                    # between specimen and cell suspension
                    else:
                        output_index = 0
                        while row.get(
                            f"{link_index_prefix}outputs.{str(output_index)}"
                        ):
                            suspension_id = row.get(
                                f"{link_index_prefix}outputs.{str(output_index)}"
                            )
                            suspension = self.cell_suspensions.get(suspension_id)

                            # Link together specimen and suspension
                            input_index = 0
                            while row.get(
                                f"{link_index_prefix}inputs.{str(input_index)}"
                            ):
                                specimen_id = row.get(
                                    f"{link_index_prefix}inputs.{str(input_index)}"
                                )
                                suspension.set_specimen_from_organism(
                                    self.specimens.get(specimen_id)
                                )

                                input_index += 1

                            output_index += 1

                if (
                    attempted_new_library.is_complete()
                    and attempted_new_library not in self.libraries
                ):
                    self.libraries.append(attempted_new_library)

                link_index += 1
