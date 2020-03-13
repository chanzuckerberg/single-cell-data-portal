import json
import os
import sys
import threading
import requests

from dcp_prototype.backend.wrangling.migrations.common.metadata import (
    MetadataBase,
    MetadataArtifact,
    MetadataProject,
    MetadataContributor,
)
from dcp_prototype.backend.wrangling.migrations.common.utils import (
    set_attribute_value,
    get_entity_type,
    append_value_to_attribute,
    append_unique_value_to_attribute,
)


class DatasetMetadata:
    """
    This class performs the mapping from DCP 1.0 metadata to the new artifact schema. All the DCP 1.0 data object are
    first added with the add_entity function. When all data is gathered, the process function is called. The data
    dictionary can then be retrieved with to_dict.
    """

    def __init__(self):
        self.artifact = MetadataArtifact()
        self.project = MetadataProject()
        self.entity_data = {}
        self.missing = {}
        self.lock = threading.Lock()

    def add_entity(self, file_source, object_data):
        """Add a dcp1 entity to the entity_data dict"""
        with self.lock:
            entity_type = get_entity_type(file_source, object_data)
            self.entity_data.setdefault(entity_type, []).append(object_data)

    def process(self):
        """Process the entity data and produce a new project for the artifact"""
        self.process_project()
        self.process_donor_organism()
        self.process_specimen_from_organism()
        self.process_sequencing_protocol()
        self.process_dissociation_protocol()
        self.process_library_preparation_protocol()
        self.process_enrichment_protocols()
        self.process_cell_suspension()
        self.process_organoid()
        self.process_cell_line()
        self.process_from_hca_server()

    def process_project(self):
        """Process the one and only project entity"""
        data = self.entity_data.get("project", [])
        if len(data) != 1:
            print("Error, only one project expected, found ", len(data))
            return False

        data = data[0]
        self.artifact.projects.append(self.project)
        mapping = (
            (("project_core", "project_short_name"), "title", set_attribute_value),
            (("project_core", "protocols_io_doi"), "protocols_io_doi", append_value_to_attribute),
            (("provenance", "document_id"), "id", set_attribute_value),
            (("publications", 0, "title"), "publication_title", set_attribute_value),
            (("publications", 0, "doi"), "publication_doi", set_attribute_value),
            (("publications", 0, "pmid"), "publication_pmid", set_attribute_value),
            (("array_express_accessions", "*"), "array_express_accessions", append_value_to_attribute),
            (("insdc_project_accessions", "*"), "insdc_project_accessions", append_value_to_attribute),
            (("geo_series_accessions", "*"), "geo_series_accessions", append_value_to_attribute),
            (("biostudies_accessions", "*"), "biostudies_accessions", append_value_to_attribute),
        )

        self.process_mappings("project", data, self.project, mapping, self.missing)

        mapping = (
            # The DCP1 has much more information about contributors
            ("name", "name", set_attribute_value),
            ("institution", "institution", set_attribute_value),
        )
        for cdata in data.get("contributors", []):
            contributor = MetadataContributor()
            self.process_mappings("project/contributors", cdata, contributor, mapping, self.missing)
            self.project.contributors.append(contributor)

    def process_donor_organism(self):
        mapping = (
            (("genus_species", "*", "ontology_label"), "donor_species", append_unique_value_to_attribute),
            (
                ("development_stage", "ontology_label"),
                "donor_development_stages_at_collection",
                append_unique_value_to_attribute,
            ),
            (("diseases", "*", "ontology_label"), "donor_diseases", append_unique_value_to_attribute),
        )
        for data in self.entity_data.get("donor_organism", []):
            self.process_mappings("donor_organism", data, self.project, mapping, self.missing)

        # This is a special case for one project that lacks an ontology_label for genus_species.
        # An alternative approach to solve this problem is to provide a list of source mappings for
        # each destination mapping, in order of preference.  Then iterate through each source mapping
        # until a match is found.  However, since this only impacts one attribute for one project,
        # it's simpler to make a special case.
        if not self.project.donor_species:
            mapping = ((("genus_species", "*", "text"), "donor_species", append_unique_value_to_attribute),)
            for data in self.entity_data.get("donor_organism", []):
                self.process_mappings("donor_organism", data, self.project, mapping, self.missing)

    def process_specimen_from_organism(self):
        # specimen_from_organism
        mapping = (
            (("organ", "ontology_label"), "organs", append_unique_value_to_attribute),
            (("organ_parts", "*", "ontology_label"), "biosample_names", append_unique_value_to_attribute),
            (("diseases", "*", "ontology_label"), "biosample_diseases", append_unique_value_to_attribute),
        )
        for data in self.entity_data.get("specimen_from_organism", []):
            self.process_mappings("specimen_from_organism", data, self.project, mapping, self.missing)

    def process_sequencing_protocol(self):
        mapping = (("paired_end", "paired_end", append_unique_value_to_attribute),)
        for data in self.entity_data.get("sequencing_protocol", []):
            self.process_mappings("sequencing_protocol", data, self.project, mapping, self.missing)

    def process_dissociation_protocol(self):
        mapping = ((("method", "ontology_label"), "cell_isolation_methods", append_unique_value_to_attribute),)
        for data in self.entity_data.get("dissociation_protocol", []):
            self.process_mappings("dissociation_protocol", data, self.project, mapping, self.missing)

    def process_library_preparation_protocol(self):
        mapping = (
            ("nucleic_acid_source", "nucleic_acid_sources", append_unique_value_to_attribute),
            (
                ("input_nucleic_acid_molecule", "ontology_label"),
                "input_nucleic_acid_molecules",
                append_unique_value_to_attribute,
            ),
            (
                ("library_construction_method", "ontology_label"),
                "library_construction_methods",
                append_unique_value_to_attribute,
            ),
        )
        for data in self.entity_data.get("library_preparation_protocol", []):
            self.process_mappings("library_preparation_protocol", data, self.project, mapping, self.missing)

    def process_cell_suspension(self):
        mapping = (
            (("selected_cell_types", "*", "ontology_label"), "selected_cell_types", append_unique_value_to_attribute),
        )
        for data in self.entity_data.get("cell_suspension", []):
            self.process_mappings("cell_suspension", data, self.project, mapping, self.missing)

    def process_enrichment_protocols(self):
        mapping = (("markers", "selected_cell_markers", append_unique_value_to_attribute),)
        for data in self.entity_data.get("enrichment_protocol", []):
            self.process_mappings("enrichment_protocol", data, self.project, mapping, self.missing)

    def process_organoid(self):
        mapping = ((("model_organ", "ontology_label"), "organs", append_unique_value_to_attribute),)
        for data in self.entity_data.get("organoid", []):
            self.process_mappings("organoid", data, self.project, mapping, self.missing)

    def process_cell_line(self):
        mapping = ((("model_organ", "ontology_label"), "organs", append_unique_value_to_attribute),)

        for data in self.entity_data.get("cell_line", []):
            self.process_mappings("cell_line", data, self.project, mapping, self.missing)

    def process_from_hca_server(self):
        resp = requests.get(f"https://service.explore.data.humancellatlas.org/repository/projects/{self.project.id}")
        if resp.status_code != 200:
            print("Error", self.project.id, self.project.title)
            raise RuntimeError(f"Unable to retrieve server metadata: {self.project.id}")

        server_data = resp.json()
        category_map = dict(specimens="primary tissue", cellLines="cell line", organoids="organoid")
        mapping = (
            (
                ("samples", "*", "sampleEntityType", "*"),
                "biosample_categories",
                lambda entity, attr, value, category_map=category_map: append_unique_value_to_attribute(
                    entity, attr, category_map.get(value, "unexpected")
                ),
            ),
        )
        self.process_mappings("server_data", server_data, self.project, mapping, self.missing)

    def to_dict(self):
        """Return the artifact data as a dict"""
        jdata = json.dumps(self.artifact, default=MetadataBase.serialize)
        return json.loads(jdata)

    def process_mapping(self, source, source_tuple, dest_object, dest_attr, mapfunc):
        """This function locates attributes from source_tuple and maps them into attributes in dest_attr.  The
        purpose of this and related functions is to make writing the mappings much more clear and concise.

        source_tuple is a tuple of strings (keys).  Each key in the tuple represents a key in the 'source' dictionary.
        The function picks off the first key, and descends into the subtree represented by that key in the source
        directory.  `process_mapping` is then called recursively on the remaining keys in source_tuple. A special key
        "*" means that the source subtree is a list, and all members of the list should be processed. As a
        convenience, source_tuple may be a single string, in which case it is treated as a tuple of size one. This is
        useful when a top level key of the source dictionary is being mapped.

        dest_attr is is the name of the attribute in the dest_object.

        mapfunc transforms the source attribute to the dest attribute.
        """

        if type(source_tuple) != tuple:
            # convert non tuples to a tuple with one item
            # this is for consistent data handling
            self.process_mapping(source, (source_tuple,), dest_object, dest_attr, mapfunc)
            return

        if len(source_tuple) > 0:
            current = source_tuple[0]
            leftover = source_tuple[1:]

            if current == "*":
                # "*" signifies that the source is an array, and that we should process all entries
                for item in source:
                    self.process_mapping(item, leftover, dest_object, dest_attr, mapfunc)

            else:
                item = source[current]
                self.process_mapping(item, leftover, dest_object, dest_attr, mapfunc)

        else:
            mapfunc(dest_object, dest_attr, source)

    def process_mappings(self, label, source_data, dest_object, mapping, missing):
        """Process a collection of mappings"""
        for source_tuple, dest_attr, mapfunc in mapping:
            source = source_data
            try:
                self.process_mapping(source, source_tuple, dest_object, dest_attr, mapfunc)
            except KeyError as e:
                if source_tuple not in missing:
                    missing[source_tuple] = (label, dest_attr, e)


def combine_projects(artifact_file, new_artifact):
    """Combine a json output from a new artifact into an existing artifact file"""
    out_project = new_artifact.get("projects")[0]
    title = out_project.get("title")
    if os.path.exists(artifact_file):
        with open(artifact_file) as json_file:
            data = json.load(json_file)
            projects = data.get("projects")
            if projects is None:
                print(f"expected 'projects' in {artifact_file}")
                sys.exit(1)

            do_append = True
            for index, project in enumerate(projects):
                if project.get("title") == title:
                    print(f"replaced project {title} in {artifact_file}")
                    projects[index] = out_project
                    do_append = False
                    break
            if do_append:
                print(f"append project {title} in {artifact_file}")
                projects.append(out_project)
            out_dict = data
    else:
        out_dict = new_artifact
        print(f"write project {title} to {artifact_file}")

    return out_dict
