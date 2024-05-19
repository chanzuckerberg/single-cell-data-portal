import hashlib
import json

from backend.wmg.data.query import DeQueryCriteria

# This dictionary will act as a storage for process data based on process_id
process_storage = {}


def generate_process_id(
    criteria1: DeQueryCriteria, criteria2: DeQueryCriteria, is_group_one: bool, de_genes: list[dict[str, str]]
):
    # Convert parameters to a JSON string and encode it to bytes
    params_str = json.dumps(
        {
            "criteria1": criteria1.dict(),
            "criteria2": criteria2.dict(),
            "is_group_one": is_group_one,
            "de_genes": de_genes,
        },
        sort_keys=True,
    )
    params_bytes = params_str.encode("utf-8")

    # Generate a SHA-256 hash of the parameters
    process_id = hashlib.sha256(params_bytes).hexdigest()

    return process_id


def store_process_data(criteria1, criteria2, is_group_one, de_genes):
    process_id = generate_process_id(criteria1, criteria2, is_group_one, de_genes)
    process_storage[process_id] = (criteria1, criteria2, is_group_one, de_genes)
    return process_id


def retrieve_process_data(process_id):
    return process_storage.get(process_id)
