import logging

from backend.geneguide.pipeline.gpt_descriptions.gpt_description_generator import (
    generate_new_gpt_descriptions,
    generate_new_seo_gpt_descriptions,
)
from backend.geneguide.pipeline.ontology_tree import get_ontology_tree_builder
from backend.geneguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.geneguide.pipeline.utils import output_json_per_key

logging.basicConfig(level=logging.INFO)


def run(*, gpt_output_directory: str, gpt_seo_output_directory: str):
    ontology_tree = get_ontology_tree_builder()
    new_gpt_descriptions, new_gpt_seo_descriptions = get_new_gpt_descriptions(ontology_tree=ontology_tree)
    output_json_per_key(new_gpt_descriptions, gpt_output_directory)
    output_json_per_key(new_gpt_seo_descriptions, gpt_seo_output_directory)


def get_new_gpt_descriptions(*, ontology_tree: OntologyTreeBuilder) -> tuple[dict[str, str], dict[str, str]]:
    new_gpt_descriptions = generate_new_gpt_descriptions(all_go_term_ids=ontology_tree.id_to_name)
    new_gpt_seo_descriptions = generate_new_seo_gpt_descriptions(new_gpt_descriptions=new_gpt_descriptions)
    return new_gpt_descriptions, new_gpt_seo_descriptions
