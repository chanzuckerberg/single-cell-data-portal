import gzip
import json
import os
from collections import defaultdict
from typing import Optional

from flask import jsonify, make_response

from backend.cellguide.common.config import CellGuideConfig
from backend.cellguide.common.constants import COMPUTATIONAL_MARKER_GENES_FOLDERNAME, MARKER_GENE_PRESENCE_FILENAME
from backend.cellguide.common.providers.s3_provider import S3Provider
from backend.cellguide.common.utils import get_object_key
from backend.common.utils.ontology_parser import ontology_parser


def _defaultdict_to_dict(d):
    if isinstance(d, defaultdict):
        # Convert the defaultdict to a dict and recursively apply this function
        return {key: _defaultdict_to_dict(value) for key, value in d.items()}
    else:
        return d


def _initialize_cellguide_marker_gene_dict():
    bucket = CellGuideConfig().bucket
    s3_provider = S3Provider()

    latest_snapshot_identifier = (
        s3_provider.download_file(bucket_name=bucket, object_key=get_object_key(object="latest_snapshot_identifier"))
        .decode("utf-8")
        .strip()
    )
    compressed_data = s3_provider.download_file(
        bucket_name=bucket,
        object_key=get_object_key(
            object=f"{latest_snapshot_identifier}/{COMPUTATIONAL_MARKER_GENES_FOLDERNAME}/{MARKER_GENE_PRESENCE_FILENAME}"
        ),
    )
    marker_gene_data = json.loads(gzip.decompress(compressed_data).decode("utf-8"))
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for gene in marker_gene_data:
        for organism in marker_gene_data[gene]:
            for tissue in marker_gene_data[gene][organism]:
                for marker in marker_gene_data[gene][organism][tissue]:
                    data[organism][tissue][marker["cell_type_id"]].append(
                        {"marker_score": marker["marker_score"], "me": marker["me"], "pc": marker["pc"], "gene": gene}
                    )

    data = _defaultdict_to_dict(data)

    for organism in data:
        for tissue in data[organism]:
            for cell_type in data[organism][tissue]:
                data[organism][tissue][cell_type].sort(key=lambda x: -x["marker_score"])

    return data


def is_id(term: str) -> bool:
    return term.startswith("CL:") or term.startswith("UBERON:") or term.startswith("NCBITaxon:")


def get(organism: str, tissue: Optional[str] = None, cell_type: Optional[str] = None):
    """
    Retrieve marker gene data for a specified organism, and optionally for a specific tissue and/or cell type.

    This function handles the retrieval of marker gene data based on the provided organism, tissue, and cell type.
    It uses ontology IDs to resolve the actual names and checks if the provided IDs exist in the marker gene data.
    If the IDs are valid, it returns the corresponding data; otherwise, it returns a 404 response.

    Parameters:
        organism (str): The ontology ID or name of the organism.
        tissue (Optional[str]): The ontology ID or name of the tissue. Default is None.
        cell_type (Optional[str]): The ontology ID or name of the cell type. Default is None.

    Returns:
        Flask Response: JSON data of the marker genes for the requested parameters or a 404 error if not found.
    """
    marker_gene_data = get_marker_gene_data()

    if is_id(organism):
        organism = ontology_parser.get_term_label(organism)
    if tissue and is_id(tissue):
        tissue = ontology_parser.get_term_label(tissue)
    if tissue and cell_type and not is_id(cell_type):
        raise ValueError("cell_type must be an ID")

    if organism not in marker_gene_data:
        return make_response(jsonify({}), 404)

    if tissue and tissue not in marker_gene_data[organism]:
        return make_response(jsonify({}), 404)

    if tissue and cell_type and cell_type not in marker_gene_data[organism][tissue]:
        return make_response(jsonify({}), 404)

    if tissue and cell_type:
        return make_response(jsonify(marker_gene_data[organism][tissue][cell_type]), 200)
    elif tissue:
        return make_response(jsonify(marker_gene_data[organism][tissue]), 200)
    else:
        return make_response(jsonify(marker_gene_data[organism]), 200)


_marker_gene_data_cache = None


def get_marker_gene_data():
    global _marker_gene_data_cache
    if _marker_gene_data_cache is None:
        if os.getenv("DEPLOYMENT_STAGE") != "test":
            # Initialize the marker gene data from the latest snapshot only if not in test mode
            _marker_gene_data_cache = _initialize_cellguide_marker_gene_dict()
        else:
            # Initialize an empty structure if in test mode
            _marker_gene_data_cache = {}
    return _marker_gene_data_cache
