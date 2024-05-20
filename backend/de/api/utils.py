import openai

from backend.common.utils.ontology_parser import ontology_parser
from backend.de.api.config import DeConfig
from backend.wmg.data.query import DeQueryCriteria


def interpret_de_results(
    criteria1: DeQueryCriteria,
    criteria2: DeQueryCriteria,
    genes1: list[str],
    genes2: list[str],
) -> str:
    prompt = _craft_de_interpretation_prompt(criteria1, criteria2, genes1, genes2)

    messages = [{"role": "user", "content": prompt}]

    response = openai.ChatCompletion.create(
        model="gpt-4o-2024-05-13",
        messages=messages,
        api_key=DeConfig().openai_api_key,
    )
    return response["choices"][0]["message"]["content"], prompt


def _craft_de_interpretation_prompt(
    criteria1: DeQueryCriteria, criteria2: DeQueryCriteria, genes1: list[str], genes2: list[str]
) -> str:
    def get_term_labels(ontology_term_ids: list[str]) -> list[str]:
        return [ontology_parser.get_term_label(term_id) for term_id in ontology_term_ids]

    def format_criteria(criteria: DeQueryCriteria) -> str:
        parts = [
            f"Organism: {ontology_parser.get_term_label(criteria.organism_ontology_term_id)}",
            (
                f"Tissues: {', '.join(get_term_labels(criteria.tissue_ontology_term_ids))}"
                if criteria.tissue_ontology_term_ids
                else ""
            ),
            (
                f"Cell Types: {', '.join(get_term_labels(criteria.cell_type_ontology_term_ids))}"
                if criteria.cell_type_ontology_term_ids
                else ""
            ),
            (
                f"Diseases: {', '.join(get_term_labels(criteria.disease_ontology_term_ids))}"
                if criteria.disease_ontology_term_ids
                else ""
            ),
            (
                f"Ethnicities: {', '.join(get_term_labels(criteria.self_reported_ethnicity_ontology_term_ids))}"
                if criteria.self_reported_ethnicity_ontology_term_ids
                else ""
            ),
            (
                f"Sexes: {', '.join(get_term_labels(criteria.sex_ontology_term_ids))}"
                if criteria.sex_ontology_term_ids
                else ""
            ),
        ]
        return ", ".join([part for part in parts if part])

    formatted_criteria1 = format_criteria(criteria1)
    formatted_criteria2 = format_criteria(criteria2)
    formatted_genes1 = "\n  ".join([f"{i+1}. {gene}" for i, gene in enumerate(genes1)])
    formatted_genes2 = "\n  ".join([f"{i+1}. {gene}" for i, gene in enumerate(genes2)])

    prompt = f"""
Please analyze the following differential expression (DE) results and provide a detailed interpretation focusing on
biologically interesting and relevant signals, such as pathway enrichment. 

- **Group 1 Criteria**: {formatted_criteria1}
- **Group 2 Criteria**: {formatted_criteria2}
- **Top Upregulated Genes for Group 1 (Downregulated for Group 2)**:
{formatted_genes1}
- **Top Downregulated Genes for Group 1 (Upregulated for Group 2)**:
{formatted_genes2}

Please include in your analysis any significant pathways, biological processes, or functional annotations that are enriched
in these top differentially expressed genes. Also, describe any notable patterns or trends that emerge from these results,
especially given the biological context of the query.
"""
    return prompt.strip()
