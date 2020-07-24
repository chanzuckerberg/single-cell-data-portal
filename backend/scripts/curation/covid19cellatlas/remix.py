import argparse
import json
import string

import scanpy as sc
import yaml

import ontology

REPLACE_SUFFIX = "_original"
ONTOLOGY_SUFFIX = "_ontology_term_id"


def is_curie(value):
    """Return True iff the value is an OBO-id CURIE like EFO:000001"""
    return ":" in value and all(c in string.digits for c in value.split(":")[1])


def is_ontology_field(field_name):
    """Return True iff the field_name is an ontology field like tissue_ontology_term_id"""
    return field_name.endswith(ONTOLOGY_SUFFIX)


def get_label_field_name(field_name):
    """Get the associated label field from an ontology field, assay_ontology_term_id --> assay"""
    return field_name[: -len(ONTOLOGY_SUFFIX)]


def safe_add_field(adata_attr, field_name, field_value):
    """Add a field and value to an AnnData, but don't clobber an exising value."""

    if isinstance(field_value, list) and field_value and isinstance(field_value[0], dict):
        field_value = json.dumps(field_value)
    if field_name in adata_attr:
        adata_attr[field_name + REPLACE_SUFFIX] = adata_attr[field_name]
    adata_attr[field_name] = field_value


def get_curie_and_label(maybe_curie):
    """Given a string that might be a curie, return a (curie, label) pair"""
    if not is_curie(maybe_curie):
        return ("", maybe_curie)
    else:
        return (maybe_curie, ontology.get_ontology_label(maybe_curie))


def remix_uns(adata, uns_config):
    """Add fields from the config to adata.uns"""
    for field_name, field_value in uns_config.items():

        if is_ontology_field(field_name):
            # If it's an ontology field, look it up
            label_field_name = get_label_field_name(field_name)
            ontology_term, ontology_label = get_curie_and_label(field_value)
            safe_add_field(adata.uns, field_name, ontology_term)
            safe_add_field(adata.uns, label_field_name, ontology_label)
        else:
            safe_add_field(adata.uns, field_name, field_value)


def remix_obs(adata, obs_config):
    """Add fields from the config to adata.obs"""

    for field_name, field_value in obs_config.items():

        if isinstance(field_value, dict):
            # If the value is a dict, that means we are supposed to map from an existing column to the new one
            source_column, column_map = next(iter(field_value.items()))

            for key in column_map:
                if key not in adata.obs[source_column].unique():
                    print(f'WARNING: key {key} not in adata.obs["{source_column}"]')

            for value in adata.obs[source_column].unique():
                if value not in column_map:
                    print(f'WARNING: Value {value} in adata.obs["{source_column}"] not in translation dict')
            if is_ontology_field(field_name):
                ontology_term_map, ontology_label_map = {}, {}
                print("\n", field_name, "\n")
                for original_value, maybe_curie in column_map.items():
                    curie, label = get_curie_and_label(maybe_curie)
                    ontology_term_map[original_value] = curie
                    ontology_label_map[original_value] = label
                    print(";".join([original_value, curie, label]))

                ontology_column = adata.obs[source_column].replace(ontology_term_map, inplace=False)
                label_column = adata.obs[source_column].replace(ontology_label_map, inplace=False)

                safe_add_field(adata.obs, field_name, ontology_column)
                safe_add_field(adata.obs, get_label_field_name(field_name), label_column)
        else:
            if is_ontology_field(field_name):
                # If it's an ontology field, look it up
                label_field_name = get_label_field_name(field_name)
                ontology_term, ontology_label = get_curie_and_label(field_value)
                safe_add_field(adata.obs, field_name, ontology_term)
                safe_add_field(adata.obs, label_field_name, ontology_label)
            else:
                safe_add_field(adata.obs, field_name, field_value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-h5ad", required=True)
    parser.add_argument("--remix-config", required=True)
    parser.add_argument("--output-filename", required=True)
    args = parser.parse_args()

    config = yaml.load(open(args.remix_config), Loader=yaml.FullLoader)
    adata = sc.read_h5ad(args.source_h5ad)
    remix_uns(adata, config["uns"])
    remix_obs(adata, config["obs"])

    adata.write_h5ad(args.output_filename, compression="gzip")


if __name__ == "__main__":
    main()
