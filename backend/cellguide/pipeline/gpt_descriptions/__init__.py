import logging

from backend.cellguide.pipeline.gpt_descriptions.gpt_description_generator import (
    generate_new_gpt_descriptions,
    upload_descriptions_to_s3,
)
from backend.cellguide.pipeline.ontology_tree import get_ontology_tree_builder
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.wmg.api.wmg_api_config import WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.snapshot import load_snapshot

logging.basicConfig(level=logging.INFO)


def run():
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    ontology_tree = get_ontology_tree_builder(snapshot=snapshot)
    new_gpt_descriptions = generate_new_gpt_descriptions(ontology_tree=ontology_tree)
    upload_descriptions_to_s3(new_gpt_descriptions=new_gpt_descriptions)


def get_new_gpt_descriptions(ontology_tree: OntologyTreeBuilder) -> dict[str, str]:
    return generate_new_gpt_descriptions(all_cell_type_ids_in_corpus=ontology_tree.all_cell_type_ids_in_corpus)
