import json
import os
import pathlib
from dataclasses import dataclass, is_dataclass

import pandas as pd

from backend.cellguide.pipeline.constants import ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME
from backend.common.utils.dataclass import convert_dataclass_to_dict


def output_json(data, path):
    """
    Output JSON while handling the case where `data` is an instance of a dataclass

    Arguments
    ---------
    data - a dict or a dataclass
    path - str
    """

    if is_dataclass(data):
        data = convert_dataclass_to_dict_and_strip_nones(data)

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f)


def convert_dataclass_to_dict_and_strip_nones(data):
    """
    A wrapper for dataclass-to-dict conversion that deletes keys with None as values
    """
    return _remove_none_values(convert_dataclass_to_dict(data))


def _remove_none_values(data):
    """
    Recursively remove keys with None values.

    Arguments
    ---------
    data - dict
    """
    if isinstance(data, dict):
        return {k: _remove_none_values(v) for k, v in data.items() if v is not None}
    elif isinstance(data, list):
        return [_remove_none_values(i) for i in data]
    else:
        return data


@dataclass
class GeneMetadata:
    gene_id_to_name: dict[str, str]
    gene_id_to_symbol: dict[str, str]


def get_gene_id_to_name_and_symbol() -> GeneMetadata:
    """
    This function reads a file containing gene metadata and returns a GeneMetadata object.
    The GeneMetadata object contains two dictionaries: gene_id_to_name and gene_id_to_symbol.
    gene_id_to_name is a dictionary where the keys are Ensembl GeneIDs and the values are gene descriptions.
    gene_id_to_symbol is a dictionary where the keys are Ensembl GeneIDs and the values are gene symbols.

    Returns
    -------
    GeneMetadata
        An object containing two dictionaries: gene_id_to_name and gene_id_to_symbol.
    """

    file_path = pathlib.Path(__file__).absolute().parent.joinpath("fixtures", ENSEMBL_GENE_ID_TO_DESCRIPTION_FILENAME)
    gene_metadata = pd.read_csv(file_path, sep="\t")
    gene_id_to_name = gene_metadata.set_index("Ensembl GeneIDs")["Description"].to_dict()
    gene_id_to_symbol = gene_metadata.set_index("Ensembl GeneIDs")["Symbols"].to_dict()
    return GeneMetadata(gene_id_to_name=gene_id_to_name, gene_id_to_symbol=gene_id_to_symbol)
