from flask import jsonify, make_response

from backend.cellguide.api.common.data import get_marker_gene_data


def get():
    """
    Retrieve all CellGuide marker gene data.

    This function handles the retrieval of all CellGuide marker gene data.

    Returns:
        Flask Response: JSON data of the marker genes.
        The response structure is a nested dictionary with the following structure:
        - organisms --> tissues --> cell types --> list(marker genes).

        Organisms and tissues are labels, cell types are IDs, and marker genes are dictionaries with the following keys:
        - `marker_score`: The score indicating the strength of the marker gene for the cell type.
        - `me`: Mean expression of the gene across the cells of the specified type.
        - `pc`: Percentage of cells within the specified type that express the gene.
        - `gene`: The gene symbol associated with the marker gene data.
    """
    return make_response(jsonify(get_marker_gene_data()), 200)
