import json

import numpy as np
import openai

from backend.de.config import DeConfig
from backend.de.data.snapshot import DeSnapshot


def get_embedding(text, model="text-embedding-ada-002"):
    text = text.replace("\n", " ")
    return np.array(
        openai.Embedding.create(input=[text], model=model, api_key=DeConfig().openai_api_key)["data"][0]["embedding"]
    )


def find_most_similar_key(key, embeddings):
    query_emb = get_embedding(key)
    corrs = np.array([np.corrcoef(query_emb, embeddings[k])[0, 1] for k in embeddings])
    return np.array(list(embeddings.keys()))[np.argmax(corrs)]


def translate_to_query_criteria(
    categories_to_values1: dict[str, list[str]], categories_to_values2: dict[str, list[str]], snapshot: DeSnapshot
):
    unique_vals_embeddings = snapshot.metadata_value_embeddings
    key_embeddings = unique_vals_embeddings["key_embeddings"]

    categories_translated = []
    values_translated = []
    for category in categories_to_values1:
        value = categories_to_values1[category]
        cat = find_most_similar_key(category, key_embeddings)
        categories_translated.append(cat)
        values_translated.append([find_most_similar_key(v, unique_vals_embeddings[cat]) for v in value])
    query_criteria1 = dict(
        zip([i + "s" if i != "organism_ontology_term_id" else i for i in categories_translated], values_translated)
    )

    categories_translated = []
    values_translated = []
    for category in categories_to_values2:
        value = categories_to_values2[category]
        cat = find_most_similar_key(category, key_embeddings)
        categories_translated.append(cat)
        values_translated.append([find_most_similar_key(v, unique_vals_embeddings[cat]) for v in value])
    query_criteria2 = dict(
        zip([i + "s" if i != "organism_ontology_term_id" else i for i in categories_translated], values_translated)
    )
    return query_criteria1, query_criteria2


def get_query_from_user_input(query: str, snapshot: DeSnapshot):
    # Step 1: send the conversation and available functions to GPT
    messages = [{"role": "user", "content": query}]
    functions = [
        {
            "name": "translate_to_query_criteria",
            "description": "Translate user query description into query criteria for differential expression.",
            "parameters": {
                "type": "object",
                "properties": {
                    "categories_to_values1": {
                        "type": "object",
                        "description": """dict[str,list[str]] - The metadata categories that will be queried mapped to the query values for population 1. 
        For example:
        {'disease': ['normal', "Alzheimer's"],
         'cell type': ['astrocyte'],
         'sex': ['male'],
         'tissue': ['brain']}  """,
                        "properties": {
                            "category": {
                                "type": "array",
                                "description": 'array of values to query for, e.g. ["normal", "Alzheimer\'s"]',
                                "items": {"type": "string"},
                            }
                        },
                    },
                    "categories_to_values2": {
                        "type": "object",
                        "description": """dict[str,list[str]] - The metadata categories that will be queried mapped to the query values for population 2. 
        For example:
        {'disease': ['normal', "Alzheimer's"],
         'cell type': ['astrocyte'],
         'sex': ['male'],
         'tissue': ['brain']}  """,
                        "properties": {
                            "category": {
                                "type": "array",
                                "description": 'array of values to query for, e.g. ["normal", "Alzheimer\'s"]',
                                "items": {"type": "string"},
                            }
                        },
                    },
                },
                "required": ["categories_to_values1", "categories_to_values2"],
            },
        }
    ]
    response = openai.ChatCompletion.create(
        model="gpt-3.5-turbo-0613",
        messages=messages,
        functions=functions,
        function_call="auto",
        api_key=DeConfig().openai_api_key,
    )
    response_message = response["choices"][0]["message"]

    # Step 2: check if GPT wanted to call a function
    if response_message.get("function_call"):
        # Step 3: call the function
        # Note: the JSON response may not always be valid; be sure to handle errors
        available_functions = {
            "translate_to_query_criteria": translate_to_query_criteria,
        }  # only one function in this example, but you can have multiple
        function_name = response_message["function_call"]["name"]
        fuction_to_call = available_functions[function_name]
        function_args = json.loads(response_message["function_call"]["arguments"])
        query1, query2 = fuction_to_call(
            function_args.get("categories_to_values1"),
            function_args.get("categories_to_values2"),
            snapshot,
        )
        query1 = json.loads(json.dumps(query1).replace('["all"]', "[]"))
        query2 = json.loads(json.dumps(query2).replace('["all"]', "[]"))
        if "organism_ontology_term_id" not in query1:
            query1["organism_ontology_term_id"] = "NCBITaxon:9606"
        else:
            query1["organism_ontology_term_id"] = query1["organism_ontology_term_id"][0]
        if "organism_ontology_term_id" not in query2:
            query2["organism_ontology_term_id"] = "NCBITaxon:9606"
        else:
            query2["organism_ontology_term_id"] = query2["organism_ontology_term_id"][0]

        return query1, query2
