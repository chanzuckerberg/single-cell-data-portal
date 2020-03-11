import json
import os
import sys
import threading
from urllib.parse import urlparse

import dcp_prototype
import jsonschema
import pkg_resources

from dcp_prototype.backend.wrangling.migrations.utils.migration_schema_classes import (
    MetadataBase,
    MetadataArtifact,
    MetadataProject,
    MetadataContributor,
)


def set_string(entity, attr, value):
    """Used in process_mapping:  set a string value to an attribute in an object"""
    setattr(entity, attr, str(value))


def unique_append_string(entity, attr, value):
    """Used in process_mapping:  append a value to a list attribute in an object, if it does not exist"""
    if value not in getattr(entity, attr):
        getattr(entity, attr).append(value)


def append_string(entity, attr, value):
    """Used in process_mapping:  append a value to a list attribute in an object"""
    getattr(entity, attr).append(value)


def process_mapping(source, source_tuple, dest_object, dest_attr, mapfunc):
    """This function locates attributes from source_tuple and maps them into attributes in dest_attr.  The purpose of
    this and related functions is to make writing the mappings much more clear and concise.

    source_tuple is a tuple of strings (keys).  Each key in the tuple represents a key in the 'source' dictionary.
    The function picks off the first key, and descends into the subtree represented by that key in the source
    directory.  `process_mapping` is then called recursively on the remaining keys in source_tuple. A special key "*"
    means that the source subtree is a list, and all members of the list should be processed. As a convenience,
    source_tuple may be a single string, in which case it is treated as a tuple of size one. This is useful when a
    top level key of the source dictionary is being mapped.

    dest_attr is is the name of the attribute in the dest_object.

    mapfunc transforms the source attribute to the dest attribute.
    """

    if type(source_tuple) != tuple:
        # convert non tuples to a tuple with one item
        # this is for consistent data handling
        process_mapping(source, (source_tuple,), dest_object, dest_attr, mapfunc)
        return

    if len(source_tuple) > 0:
        current = source_tuple[0]
        leftover = source_tuple[1:]

        if current == "*":
            # "*" signifies that the source is an array, and that we should
            # process all entries
            for item in source:
                process_mapping(item, leftover, dest_object, dest_attr, mapfunc)

        else:
            item = source[current]
            process_mapping(item, leftover, dest_object, dest_attr, mapfunc)

    else:
        mapfunc(dest_object, dest_attr, source)


def process_mappings(label, source_data, dest_object, mapping, missing):
    """Process a collection of mappings"""
    for source_tuple, dest_attr, mapfunc in mapping:
        source = source_data
        try:
            process_mapping(source, source_tuple, dest_object, dest_attr, mapfunc)
        except KeyError as e:
            if source_tuple not in missing:
                missing[source_tuple] = (label, dest_attr, e)


class DatasetMetadata(object):
    """This class performs the mapping from dcp1 to the new artifact schema.
    All the dcp1 data object are first added with the add_entity function.
    When all data is gathered, the process function is called.
    The data dictionary can then be retrieved with to_dict.
    """

    def __init__(self):
        self.artifact = MetadataArtifact()
        self.project = MetadataProject()
        self.entity_data = {}
        self.missing = {}
        self.lock = threading.Lock()

    def validate(self, data):
        """Validate the data with the current artifact schema"""
        try:
            schema_file = pkg_resources.resource_filename(dcp_prototype.__name__, "backend/v0.0.0.json")
            schema_data = open(schema_file).read()
            schema = json.loads(schema_data)
        except FileNotFoundError:
            # FIXME, the schema is not yet committed to main.
            # Validate should fail if the data cannot be validated.
            print(f"Warning, the schema could not be validated: {schema_file}")
            return False

        try:
            jsonschema.validate(data, schema)
            return True
        except jsonschema.ValidationError as e:
            print("validation error: ", e)
            return False

    def add_entity(self, file_source, object_data):
        """Add a dcp1 entity to the entity_data dict"""
        with self.lock:
            entity_type = get_entity_type(file_source, object_data)
            self.entity_data.setdefault(entity_type, []).append(object_data)

    def process_project(self):
        """Process the one and only project entity"""
        data = self.entity_data.get("project", [])
        if len(data) != 1:
            print("Error, only one project expected, found ", len(data))
            return False

        data = data[0]
        self.artifact.projects.append(self.project)
        mapping = (
            (("project_core", "project_short_name"), "title", set_string),
            (("project_core", "protocols_io_doi"), "protocols_io_doi", append_string),
            (("provenance", "document_id"), "id", set_string),
            (("publications", 0, "title"), "publication_title", set_string),
            (("publications", 0, "doi"), "publication_doi", set_string),
            (("publications", 0, "pmid"), "publication_pmid", set_string),
            (("array_express_accessions", "*"), "array_express_accessions", append_string),
            (("insdc_project_accessions", "*"), "insdc_project_accessions", append_string),
            (("geo_series_accessions", "*"), "geo_series_accessions", append_string),
            (("biostudies_accessions", "*"), "biostudies_accessions", append_string),
        )

        process_mappings("project", data, self.project, mapping, self.missing)

        mapping = (
            # The DCP1 has much more information about contributors
            ("name", "name", set_string),
            ("institution", "institution", set_string),
        )
        for cdata in data.get("contributors", []):
            contributor = MetadataContributor()
            process_mappings("project/contributors", cdata, contributor, mapping, self.missing)
            self.project.contributors.append(contributor)

    def process_donor_organism(self):
        mapping = (
            (("genus_species", "*", "ontology_label"), "donor_species", unique_append_string),
            (("development_stage", "ontology_label"), "donor_development_stages_at_collection", unique_append_string),
            (("diseases", "*", "ontology_label"), "donor_diseases", unique_append_string),
        )
        for data in self.entity_data.get("donor_organism", []):
            process_mappings("donor_organism", data, self.project, mapping, self.missing)

    def process_specimen_from_organism(self):
        # specimen_from_organism
        mapping = (
            (("biomaterial_core", "biomaterial_name"), "biosample_names", unique_append_string),
            (("organ", "ontology_label"), "organs", unique_append_string),
            (("organ_parts", "*", "ontology_label"), "biosample_categories", unique_append_string),
            (("diseases", "*", "ontology_label"), "biosample_diseases", unique_append_string),
        )
        for data in self.entity_data.get("specimen_from_organism", []):
            process_mappings("specimen_from_organism", data, self.project, mapping, self.missing)

    def process_sequencing_protocol(self):
        mapping = (("paired_end", "paired_end", unique_append_string),)
        for data in self.entity_data.get("sequencing_protocol", []):
            process_mappings("sequencing_protocol", data, self.project, mapping, self.missing)

    def process_dissociation_protocol(self):
        mapping = ((("method", "ontology_label"), "cell_isolation_methods", unique_append_string),)
        for data in self.entity_data.get("dissociation_protocol", []):
            process_mappings("dissociation_protocol", data, self.project, mapping, self.missing)

    def process_library_preparation_protocol(self):
        mapping = (
            ("nucleic_acid_source", "nucleic_acid_sources", unique_append_string),
            (("input_nucleic_acid_molecule", "ontology_label"), "input_nucleic_acid_molecules", unique_append_string),
            (("library_construction_method", "ontology_label"), "library_construction_methods", unique_append_string),
        )
        for data in self.entity_data.get("library_preparation_protocol", []):
            process_mappings("library_preparation_protocol", data, self.project, mapping, self.missing)

    def process_cell_suspension(self):
        mapping = ((("selected_cell_types", "*", "ontology_label"), "selected_cell_types", unique_append_string),)
        for data in self.entity_data.get("cell_suspension", []):
            process_mappings("cell_suspension", data, self.project, mapping, self.missing)

    def process_enrichment_protocols(self):
        mapping = (("markers", "selected_cell_markers", unique_append_string),)
        for data in self.entity_data.get("cell_suspension", []):
            process_mappings("cell_suspension", data, self.project, mapping, self.missing)

    def process(self):
        """Process the entity data and produce a new project for the artifact"""
        self.process_project()
        self.process_donor_organism()
        self.process_specimen_from_organism()
        self.process_sequencing_protocol()
        self.process_dissociation_protocol()
        self.process_library_preparation_protocol()
        self.process_enrichment_protocols()

    def to_dict(self):
        """Return the artifact data as a dict"""
        jdata = json.dumps(self.artifact, default=MetadataBase.serialize)
        return json.loads(jdata)


def get_entity_type(object_name, object_data):
    """Return the entity_type from the object_data"""
    desc = object_data.get("describedBy")
    if desc:
        path = urlparse(desc).path
        entity_type = os.path.basename(path)
        return entity_type
    else:
        raise RuntimeError(f"object {object_name} did not have a describedBy field")


def combine_projects(artifact_file, new_artifact):
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
