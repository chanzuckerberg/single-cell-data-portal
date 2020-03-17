"""
Simple script to convert a DCP/1 matrix service loom to the DCP/2 cell
metadata.

Modifies the loom file in place and ends up dropping quite a few of the
fields from DCP/1.
"""

import argparse
import loompy

# Map from dcp/1 matrix metadata fields to dcp/2
METADATA_MAP = {
    "CellID": "CellID",
    "barcode": None,
    "bundle_uuid": None,
    "bundle_version": None,
    "cell_suspension.provenance.document_id": None,
    "cell_suspension.genus_species.ontology": "biosample_preparation.donor_species",
    "cell_suspension.genus_species.ontology_label": None,
    "derived_organ_label": None,
    "derived_organ_ontology": "biosample_preparation.organ",
    "derived_organ_parts_label": None,
    "derived_organ_parts_ontology": None,
    "donor_organism.development_stage.ontology": "biosample_preparation.donor_development_stage_at_collection",
    "donor_organism.development_stage.ontology_label": None,
    "donor_organism.diseases.ontology": "biosample_preparation.donor_diseases",
    "donor_organism.diseases.ontology_label": None,
    "donor_organism.human_specific.ethnicity.ontology": "biosample_preparation.donor_ethnicity",
    "donor_organism.human_specific.ethnicity.ontology_label": None,
    "donor_organism.is_living": None,
    "donor_organism.provenance.document_id": "biosample_preparation.donor_accession",
    "donor_organism.sex": "biosample_preparation.donor_sex",
    "emptydrops_is_cell": None,
    "dss_bundle_fqid": None,
    "file_uuid": None,
    "file_version": None,
    "genes_detected": None,
    "library_preparation_protocol.end_bias": None,
    "library_preparation_protocol.input_nucleic_acid_molecule.ontology": (
        "library_preparation_protocol.nucleic_acid_source"
    ),
    "library_preparation_protocol.input_nucleic_acid_molecule.ontology_label": None,
    "library_preparation_protocol.library_construction_method.ontology": (
        "library_preparation_protocol.library_construction_method"
    ),
    "library_preparation_protocol.library_construction_method.ontology_label": None,
    "library_preparation_protocol.provenance.document_id": None,
    "library_preparation_protocol.strand": None,
    "project.project_core.project_short_name": None,
    "specimen_from_organism.organ.ontology": None,
    "specimen_from_organism.organ.ontology_label": None,
    "specimen_from_organism.organ_parts.ontology": None,
    "specimen_from_organism.organ_parts.ontology_label": None,
    "specimen_from_organism.provenance.document_id": None,
    "total_umis": None,
    "project.project_core.project_title": "project.project_title",
    "project.provenance.document_id": None,
    "analysis_protocol.protocol_core.protocol_id": None,
    "analysis_protocol.provenance.document_id": None,
    "analysis_working_group_approval_status": None,
}


def get_category(org):
    """Get the "biosample category" based on the organ. The organ field in
    DCP/1 matrix metadata contained information about the biosample category.
    """

    if "cell line" in org:
        return "cell line"
    elif "organoid" in org:
        return "organoid"
    return "primary"


def strip_category(org):
    """Remove the category annotation stuck to the end of organ fields
    in dcp/1 metadata.
    """
    return org.rstrip(" (organoid)").rstrip(" (cell line)")


def convert_loom_metadata(loom_path):
    """Convert the metadata in loom_path from dcp/1 to dcp/2"""

    ds = loompy.connect(loom_path)

    for dcp1, dcp2 in METADATA_MAP.items():
        if dcp1 not in ds.ca:
            continue
        if not dcp2:
            continue
        ds.ca[dcp2] = ds.ca[dcp1]

    for dcp1 in METADATA_MAP:
        if dcp1 not in ds.ca:
            continue
        if dcp1 != METADATA_MAP[dcp1]:
            delattr(ds.ca, dcp1)

    categories = [get_category(o) for o in ds.ca["biosample_preparation.organ"]]
    ds.ca["biosample_preparation.biosample_category"] = categories

    organs = [strip_category(o) for o in ds.ca["biosample_preparation.organ"]]
    ds.ca["biosample_preparation.organ"] = organs

    ds.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("loom_file", help="Loom file to modify in place to have dcp/2 metadata")
    args = parser.parse_args()
    convert_loom_metadata(args.loom_file)
