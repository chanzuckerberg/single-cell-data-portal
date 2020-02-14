from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_project import (  # noqa
    OldProject,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_sequencing_protocol import (  # noqa
    OldSequencingProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_library_prep_protocol import (  # noqa
    OldLibraryPrepProtocol,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_cell_suspension import (  # noqa
    OldCellSuspension,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_sequence_file import (  # noqa
    OldSequenceFile,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_donor_organism import (  # noqa
    OldDonorOrganism,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_specimen_from_organism import (  # noqa
    OldSpecimenFromOrganism,
)
from dcp_prototype.backend.wrangling.migrations.metadata_schema_representation.old_entities.old_analysis_file import (  # noqa
    OldAnalysisFile,
)

from dcp_prototype.backend.wrangling.migrations.utils.util import merge_dictionary_into
from pandas import DataFrame, ExcelWriter
import random
from dcp_prototype.backend.ledger.code.common.ledger_orm import (
    DBSessionMaker,
    AlignmentProtocol,
    Library,
    ProjectContributorJoin,
    SequenceFileAlignmentProtocolAnalysisFileProcessJoin,
    BiosamplePrepLibraryLibraryPrepProtocolProcessJoin,
    LibrarySequenceFileSequencingProtocolProcessJoin,
)
from dcp_prototype.backend.wrangling.migrations.utils.constants import (
    SS2_ALIGNMENT_PROTOCOL,
)
from dcp_prototype.backend.wrangling.migrations.utils.id_generator import (
    hca_accession_generator,
    hca_accession_transformer,
)
import threading

ENTITY_TYPES = [
    "donor_organism",
    "specimen_from_organism",
    "cell_suspension",
    "library_preparation_protocol",
    "project",
    "sequence_file",
    "links",
    "sequencing_protocol",
    "analysis_file",
]

class OldDatasetMetadata:
    def __init__(self, sequencing_technology="ss2", s3_uri=None):
        self.sequencing_technology = sequencing_technology
        self.s3_uri = s3_uri

        self.donor_organisms = {}
        self.specimens = {}
        self.cell_suspensions = {}
        self.library_preps = {}
        self.projects = {}
        self.contributors = {}
        self.sequence_files = {}
        self.sequencing_protocols = {}
        self.analysis_files = {}

        self.libraries = []
        self.project_contributor_links = []
        self.sequence_file_analysis_file_links = []
        self.biosample_library_prep_project_links = []
        self.project_sequencing_protocol_sequence_file_links = []

        self.publish_mode = False
        self.old_metadata_structure = {}
        self.new_metadata_structure = {}

        self.project_name = ""

        self.locks = {}
        for etype in ENTITY_TYPES:
            self.locks[etype] = threading.Lock()


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

    def export_to_database(self, db_session_maker: DBSessionMaker):
        session = db_session_maker.session()

        # STEP 1: POPULATE ENTITIES

        # BiosamplePrep table population
        biosample_preps = {}
        for id, old_biosample_prep in self.cell_suspensions.items():
            biosample_prep = old_biosample_prep.convert_to_new_entity()
            biosample_preps[id] = biosample_prep
            session.add(biosample_prep)

        # LibraryPrepProtocol table population
        library_prep_protocols = {}
        for id, old_library_prep in self.library_preps.items():
            library_prep = old_library_prep.convert_to_new_entity()
            library_prep_protocols[id] = library_prep
            session.add(library_prep)

        # Project table population
        projects = {}
        for id, old_project in self.projects.items():
            project = old_project.convert_to_new_entity()
            projects[id] = project
            session.add(project)

        # Contributor table population
        contributors = {}
        for id, old_contributor in self.contributors.items():
            contributor = old_contributor.convert_to_new_entity()
            contributors[id] = contributor
            session.add(contributor)

        # SequenceProtocol table population
        sequencing_protocols = {}
        for id, old_sequencing_protocol in self.sequencing_protocols.items():
            sequencing_protocol = old_sequencing_protocol.convert_to_new_entity()
            sequencing_protocols[id] = sequencing_protocol
            session.add(sequencing_protocol)

        # SequenceFile table population
        sequence_files = {}
        for id, old_sequence_file in self.sequence_files.items():
            sequence_file = old_sequence_file.convert_to_new_entity()
            sequence_files[id] = sequence_file
            session.add(sequence_file)

        # AnalysisFile table population
        analysis_files = {}
        for id, old_analysis_file in self.analysis_files.items():
            analysis_file = old_analysis_file.convert_to_new_entity()
            analysis_files[id] = analysis_file
            session.add(analysis_file)

        # AlignmentProtocol table population
        alignment_protocols = {}
        alignment_protocol = SS2_ALIGNMENT_PROTOCOL
        alignment_protocol.id = hca_accession_generator(AlignmentProtocol.__name__)
        alignment_protocols[alignment_protocol.id] = alignment_protocol
        session.add(alignment_protocol)

        # STEP 2: POPULATE LINKS

        # Library table population
        libraries = {}
        for id, project in projects.items():
            library_id = hca_accession_transformer(Library.__name__, id)
            library = Library(id=library_id, project=project)
            libraries[id] = library
            session.add(library)

        # ProjectContributorJoin table population
        for link in self.project_contributor_links:
            project = projects.get(link[0].corresponding_old_id)
            contributor = contributors.get(link[1].name)
            id = hca_accession_generator(ProjectContributorJoin.__name__)
            project_contributor = ProjectContributorJoin(
                id=id, project=project, contributor=contributor
            )
            session.add(project_contributor)

        # SequenceFileAlignmentProtocolAnalysisFileProcessJoin table population
        skipped_count = 0
        for link in self.sequence_file_analysis_file_links:
            if link[0] is None or link[1] is None:
                skipped_count += 1
                continue
            sequence_file = sequence_files.get(link[0].corresponding_old_id)
            analysis_file = analysis_files.get(link[1].corresponding_old_id)
            alignment_protocol = random.choice(list(alignment_protocols.values()))
            id = hca_accession_generator(
                SequenceFileAlignmentProtocolAnalysisFileProcessJoin.__name__
            )
            join_object = SequenceFileAlignmentProtocolAnalysisFileProcessJoin(
                id=id,
                sequence_file=sequence_file,
                analysis_file=analysis_file,
                alignment_protocol=alignment_protocol,
            )
            session.add(join_object)
        print(
            f"Had to skip {skipped_count} objects for "
            f"SequenceFileAlignmentProtocolAnalysisFileProcessJoin."
        )

        # BiosamplePrepLibraryLibraryPrepProtocolProcessJoin table population
        for link in self.biosample_library_prep_project_links:
            biosample = biosample_preps.get(link[0].corresponding_old_id)
            library_prep = library_prep_protocols.get(link[1].corresponding_old_id)
            library = libraries.get(link[2].corresponding_old_id)
            id = hca_accession_generator(
                BiosamplePrepLibraryLibraryPrepProtocolProcessJoin.__name__
            )
            join_object = BiosamplePrepLibraryLibraryPrepProtocolProcessJoin(
                id=id,
                biosample_prep=biosample,
                library=library,
                library_prep_protocol=library_prep,
            )
            session.add(join_object)

        # LibrarySequenceFileSequencingProtocolProcessJoin table population
        skipped_count = 0
        for link in self.project_sequencing_protocol_sequence_file_links:
            if link[0] is None or link[1] is None or link[2] is None:
                skipped_count += 1
                continue
            library = libraries.get(link[0].corresponding_old_id)
            sequencing_protocol = sequencing_protocols.get(link[1].corresponding_old_id)
            sequence_file = sequence_files.get(link[2].corresponding_old_id)
            id = hca_accession_generator(
                LibrarySequenceFileSequencingProtocolProcessJoin.__name__
            )
            join_object = LibrarySequenceFileSequencingProtocolProcessJoin(
                id=id,
                library=library,
                sequence_file=sequence_file,
                sequencing_protocol=sequencing_protocol,
            )
            session.add(join_object)

        print(f"Had to skip {skipped_count} objects.")
        print(f"Beginning commit to database")
        session.commit()

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

        self.old_metadata_structure[
            "Donor Organism"
        ] = self.donor_organisms.values()
        self.old_metadata_structure[
            "Specimen From Organism"
        ] = self.specimens.values()
        self.old_metadata_structure[
            "Cell Suspensions"
        ] = self.cell_suspensions.values()
        self.old_metadata_structure[
            "Library Prep Protocols"
        ] = self.library_preps.values()
        self.old_metadata_structure["Libraries"] = self.libraries
        self.old_metadata_structure["Projects"] = self.projects.values()
        self.old_metadata_structure["Contributors"] = self.contributors.values()
        self.old_metadata_structure[
            "Sequencing Protocols"
        ] = self.sequencing_protocols.values()
        self.old_metadata_structure["Sequence Files"] = self.sequence_files.values()

        if self.old_metadata_structure["Projects"]:
            self.project_name = list(self.old_metadata_structure["Projects"])[
                0
            ].project_short_name
        else:
            self.project_name = "For some reason this project got no name."

        self.publish_mode = True

    def parse_flattened_row_of_json(self, row, entity_type):
        if entity_type == "donor_organism":
            donor_organism = OldDonorOrganism()
            donor_organism.populate_from_dcp_one_json_data_frame(row)
            with self.locks[entity_type]:
                if donor_organism.corresponding_old_id not in self.donor_organisms.keys():
                    self.donor_organisms[
                        donor_organism.corresponding_old_id
                    ] = donor_organism

        if entity_type == "specimen_from_organism":
            specimen_from_organism = OldSpecimenFromOrganism()
            specimen_from_organism.populate_from_dcp_one_json_data_frame(row)
            with self.locks[entity_type]:
                if specimen_from_organism.corresponding_old_id not in self.specimens.keys():
                    self.specimens[
                        specimen_from_organism.corresponding_old_id
                    ] = specimen_from_organism

        if entity_type == "cell_suspension":
            cell_suspension = OldCellSuspension()
            cell_suspension.populate_from_dcp_one_json_data_frame(row)
            with self.locks[entity_type]:
                if cell_suspension.corresponding_old_id not in self.cell_suspensions.keys():
                    self.cell_suspensions[
                        cell_suspension.corresponding_old_id
                    ] = cell_suspension

        if entity_type == "library_preparation_protocol":
            library_prep = OldLibraryPrepProtocol()
            library_prep.populate_from_dcp_one_json_data_frame(row)
            with self.locks[entity_type]:
                if library_prep.corresponding_old_id not in self.library_preps.keys():
                    self.library_preps[library_prep.corresponding_old_id] = library_prep

        if entity_type == "project":
            project = OldProject()
            contributors = project.populate_from_dcp_one_json_data_frame(row)
            with self.locks[entity_type]:
                if project.corresponding_old_id not in self.projects.keys():
                    self.projects[project.corresponding_old_id] = project

                    for contributor in contributors:
                        if contributor.name not in self.contributors.keys():
                            self.contributors[contributor.name] = contributor
                            self.project_contributor_links.append((project, contributor))

        if entity_type == "sequencing_protocol":
            protocol = OldSequencingProtocol()
            protocol.populate_from_dcp_one_json_data_frame(row)
            with self.locks[entity_type]:
                if protocol.corresponding_old_id not in self.sequencing_protocols.keys():
                    self.sequencing_protocols[protocol.corresponding_old_id] = protocol

        if entity_type == "sequence_file":
            sequence_file = OldSequenceFile()
            sequence_file.populate_from_dcp_one_json_data_frame(row)
            sequence_file.set_s3_uri(self.s3_uri)
            with self.locks[entity_type]:
                if sequence_file.corresponding_old_id not in self.sequence_files.keys():
                    self.sequence_files[sequence_file.corresponding_old_id] = sequence_file

        if entity_type == "analysis_file":
            analysis_file = OldAnalysisFile()
            analysis_file.populate_from_dcp_one_json_data_frame(row)
            analysis_file.set_s3_uri(self.s3_uri)
            with self.locks[entity_type]:
                if analysis_file.corresponding_old_id not in self.analysis_files.keys():
                    self.analysis_files[analysis_file.corresponding_old_id] = analysis_file

        if entity_type == "links":
            with self.locks[entity_type]:
                self.parse_links_dot_json(row)

    def parse_links_dot_json(self, row):
        link_index = 0

        while row.get(f"links.{str(link_index)}.process"):
            link_index_prefix = f"links.{str(link_index)}."

            # Link together Sequence Files and Analysis Files
            if row.get(f"{link_index_prefix}input_type") == "file":
                input_index = 0
                while row.get(f"{link_index_prefix}inputs.{str(input_index)}"):
                    seq_file_id = row.get(
                        f"{link_index_prefix}inputs.{str(input_index)}"
                    )

                    output_index = 0
                    while row.get(f"{link_index_prefix}outputs.{str(output_index)}"):
                        anal_file_id = row.get(
                            f"{link_index_prefix}outputs.{str(output_index)}"
                        )

                        self.sequence_file_analysis_file_links.append(
                            (
                                self.sequence_files.get(seq_file_id),
                                self.analysis_files.get(anal_file_id),
                            )
                        )

                        output_index += 1

                    input_index += 1

            elif (
                row.get(f"{link_index_prefix}input_type") == "biomaterial"
                and row.get(f"{link_index_prefix}output_type") == "file"
            ):
                protocol_index = 0
                while row.get(
                    f"{link_index_prefix}protocols.{str(protocol_index)}.protocol_type"
                ):
                    protocol_type = row.get(
                        f"{link_index_prefix}protocols.{str(protocol_index)}."
                        f"protocol_type"
                    )
                    protocol_id = row.get(
                        f"{link_index_prefix}protocols.{str(protocol_index)}."
                        f"protocol_id"
                    )

                    # Link together sequence files and sequence protocol and project
                    if protocol_type == "sequencing_protocol":
                        file_index = 0
                        while row.get(f"{link_index_prefix}outputs.{str(file_index)}"):
                            sequence_file_id = row.get(
                                f"{link_index_prefix}outputs.{str(file_index)}"
                            )

                            self.project_sequencing_protocol_sequence_file_links.append(
                                (
                                    random.choice(list(self.projects.values())),
                                    self.sequencing_protocols.get(protocol_id),
                                    self.sequence_files.get(sequence_file_id),
                                )
                            )

                            file_index += 1

                    # Link together biosample, library prep, and project
                    elif protocol_type == "library_preparation_protocol":
                        input_index = 0
                        while row.get(f"{link_index_prefix}inputs.{str(input_index)}"):
                            suspension_id = row.get(
                                f"{link_index_prefix}inputs.{str(input_index)}"
                            )

                            self.biosample_library_prep_project_links.append(
                                (
                                    self.cell_suspensions.get(suspension_id),
                                    self.library_preps.get(protocol_id),
                                    random.choice(list(self.projects.values())),
                                )
                            )

                            input_index += 1

                    protocol_index += 1

            # Link together donor organism -> specimen from organism -> cell suspension
            # that will form a BiosamplePrep.
            elif (
                row.get(f"{link_index_prefix}input_type") == "biomaterial"
                and row.get(f"{link_index_prefix}output_type") == "biomaterial"
            ):
                # If a link has no protocols associated, then the link is between
                # donor organism and specimen
                if row.get(
                    f"{link_index_prefix}protocols.0.protocol_type"
                ) is None or not row.get(
                    f"{link_index_prefix}protocols.0.protocol_type"
                ):
                    output_index = 0
                    while row.get(f"{link_index_prefix}outputs.{str(output_index)}"):
                        specimen_id = row.get(
                            f"{link_index_prefix}outputs.{str(output_index)}"
                        )
                        specimen = self.specimens.get(specimen_id)

                        # Link together specimen and donor organism
                        input_index = 0
                        while row.get(f"{link_index_prefix}inputs.{str(input_index)}"):
                            donor_id = row.get(
                                f"{link_index_prefix}inputs.{str(input_index)}"
                            )
                            specimen.set_donor_organism(
                                self.donor_organisms.get(donor_id)
                            )

                            input_index += 1

                        output_index += 1

                # If a link has protocols associated with it, then the link is between
                # specimen and cell suspension
                else:
                    output_index = 0
                    while row.get(f"{link_index_prefix}outputs.{str(output_index)}"):
                        suspension_id = row.get(
                            f"{link_index_prefix}outputs.{str(output_index)}"
                        )
                        suspension = self.cell_suspensions.get(suspension_id)

                        # Link together specimen and suspension
                        input_index = 0
                        while row.get(f"{link_index_prefix}inputs.{str(input_index)}"):
                            specimen_id = row.get(
                                f"{link_index_prefix}inputs.{str(input_index)}"
                            )
                            suspension.set_specimen_from_organism(
                                self.specimens.get(specimen_id)
                            )

                            input_index += 1

                        output_index += 1

            link_index += 1
