import logging

from backend.cellguide.pipeline.gpt_descriptions.gpt_description_generator import (
    generate_new_gpt_descriptions,
    generate_new_seo_gpt_descriptions,
)
from backend.cellguide.pipeline.ontology_tree import get_ontology_tree_builder
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.utils import output_json_per_key
from backend.wmg.api.wmg_api_config import WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.snapshot import load_snapshot

logging.basicConfig(level=logging.INFO)


def run(*, gpt_output_directory: str, gpt_seo_output_directory: str):
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    ontology_tree = get_ontology_tree_builder(snapshot=snapshot)
    new_gpt_descriptions, new_gpt_seo_descriptions = get_new_gpt_descriptions(ontology_tree=ontology_tree)
    output_json_per_key(new_gpt_descriptions, gpt_output_directory)
    output_json_per_key(new_gpt_seo_descriptions, gpt_seo_output_directory)


def get_new_gpt_descriptions(*, ontology_tree: OntologyTreeBuilder) -> tuple[dict[str, str], dict[str, str]]:
    new_gpt_descriptions = generate_new_gpt_descriptions(
        all_cell_type_ids_to_labels_in_corpus=ontology_tree.all_cell_type_ids_to_labels_in_corpus
    )
    new_gpt_seo_descriptions = generate_new_seo_gpt_descriptions(new_gpt_descriptions=new_gpt_descriptions)
    return new_gpt_descriptions, new_gpt_seo_descriptions
