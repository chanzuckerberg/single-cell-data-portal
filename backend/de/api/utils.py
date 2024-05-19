import openai

from backend.de.api.config import DeConfig
from backend.wmg.data.query import DeQueryCriteria


def interpret_de_results(
    criteria1: DeQueryCriteria,
    criteria2: DeQueryCriteria,
    is_group_one: bool,
    genes: list[dict[str, str]],
) -> str:
    prompt = _craft_de_interpretation_prompt(criteria1, criteria2, "1" if is_group_one else "2", genes)

    messages = [{"role": "user", "content": prompt}]

    response = openai.ChatCompletion.create(
        model="gpt-4o-2024-05-13",
        messages=messages,
        api_key=DeConfig().openai_api_key,
    )
    return response["choices"][0]["message"]["content"]


def _craft_de_interpretation_prompt(
    criteria1: DeQueryCriteria, criteria2: DeQueryCriteria, selected_group: int, genes: list[dict[str, str]]
) -> str:
    def format_criteria(criteria: DeQueryCriteria) -> str:
        parts = [
            f"Organism: {criteria.organism_ontology_term_id}",
            f"Tissues: {', '.join(criteria.tissue_ontology_term_ids)}" if criteria.tissue_ontology_term_ids else "",
            (
                f"Cell Types: {', '.join(criteria.cell_type_ontology_term_ids)}"
                if criteria.cell_type_ontology_term_ids
                else ""
            ),
            f"Diseases: {', '.join(criteria.disease_ontology_term_ids)}" if criteria.disease_ontology_term_ids else "",
            (
                f"Ethnicities: {', '.join(criteria.self_reported_ethnicity_ontology_term_ids)}"
                if criteria.self_reported_ethnicity_ontology_term_ids
                else ""
            ),
            f"Sexes: {', '.join(criteria.sex_ontology_term_ids)}" if criteria.sex_ontology_term_ids else "",
        ]
        return ", ".join([part for part in parts if part])

    formatted_criteria1 = format_criteria(criteria1)
    formatted_criteria2 = format_criteria(criteria2)
    formatted_genes = "\n  ".join(
        [
            f"{i+1}. {gene['gene_symbol']} (p-value: {gene['p_value']}, effect size: {gene['effect_size']})"
            for i, gene in enumerate(genes)
        ]
    )

    prompt = f"""
Please analyze the following differential expression (DE) results and provide a detailed interpretation focusing on biologically interesting and relevant signals, such as pathway enrichment. Factor the reported p-value and effect sizes into account. Typically, p-values < 0.001 and effect sizes > 1 are meaningful.

- **Group 1 Criteria**: {formatted_criteria1}
- **Group 2 Criteria**: {formatted_criteria2}
- **Top DE Genes for Group {selected_group}**:
  {formatted_genes}

Please include in your analysis any significant pathways, biological processes, or functional annotations that are enriched in these top differentially expressed genes. Also, describe any notable patterns or trends that emerge from these results.
"""
    return prompt.strip()
