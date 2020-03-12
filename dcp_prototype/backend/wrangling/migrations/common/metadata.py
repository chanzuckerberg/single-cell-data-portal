"""This module contains simple python classes that correspond to the artifact schema"""

import json


class MetadataBase(object):
    """Base class for the metadata objects"""

    @staticmethod
    def serialize(obj):
        """A method used by json.dumps to serialize a python object"""
        if isinstance(obj, MetadataBase):
            return obj.__dict__


class MetadataProject(MetadataBase):
    """Represents the project object defined in the Artifact Schema"""

    def __init__(self):
        """The attributes in this class must contain exactly the attributes describe in the artifact schema.
        The comment ahead of each group of attributes defines the source_tuple
        file from which that data is retrieved."""

        # project
        self.title = ""
        self.id = ""
        self.publication_title = ""
        self.publication_doi = ""
        self.publication_pmid = ""
        self.array_express_accessions = []
        self.insdc_project_accessions = []
        self.geo_series_accessions = []
        self.biostudies_accessions = []
        self.protocols_io_doi = []
        self.contributors = []

        # donor_organism
        self.donor_development_stages_at_collection = []
        self.donor_species = []
        self.donor_diseases = []

        # specimen_from_organism
        self.biosample_names = []
        self.biosample_categories = []
        self.biosample_diseases = []
        self.organs = []

        # sequencing protocol
        self.paired_end = []

        # dissociation protocol
        self.cell_isolation_methods = []

        # the number of cell_suspension files
        self.selected_cell_types = []

        # library_preparation_protocol
        self.nucleic_acid_sources = []
        self.input_nucleic_acid_molecules = []
        self.library_construction_methods = []

        # enrichment_protocol
        self.selected_cell_markers = []

        # TODO: NOT PROCESSED YET
        self.systems = []
        self.assay_categories = []
        self.cell_count = 0


class MetadataContributor(MetadataBase):
    """Defines the contributor data, as defined in the artifact schema"""

    def __init__(self):
        self.name = ""
        self.institution = ""


class MetadataArtifact(MetadataBase):
    """Defines the top level artifact data as defined by the artifact schema."""

    def __init__(self):
        self.projects = []
        self.files = []
