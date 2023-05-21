import itertools
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Tuple

import connexion
import numpy as np
import pandas as pd
from flask import jsonify
from pandas import DataFrame
from server_timing import Timing as ServerTiming

from backend.wmg.data.calculate_markers import _run_ttest
from backend.wmg.data.ontology_labels import gene_term_label, ontology_term_label
from backend.wmg.data.query import (
    FmgQueryCriteria,
    MarkerGeneQueryCriteria,
    WmgFiltersQueryCriteria,
    WmgQuery,
    WmgQueryCriteria,
    retrieve_top_n_markers,
)
from backend.wmg.data.rollup import rollup_across_cell_type_descendants
from backend.wmg.data.schemas.cube_schema import expression_summary_non_indexed_dims
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot
from backend.wmg.data.utils import depluralize, find_all_dim_option_values, find_dim_option_values, to_dict

# TODO: add cache directives: no-cache (i.e. revalidate); impl etag
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2132

DEFAULT_GROUP_BY_TERMS = ["tissue_ontology_term_id", "cell_type_ontology_term_id"]


def primary_filter_dimensions():
    snapshot: WmgSnapshot = load_snapshot()
    return jsonify(snapshot.primary_filter_dimensions)


def query():
    request = connexion.request.json
    is_rollup = request.get("is_rollup", True)
    compare = request.get("compare", None)

    if compare:
        compare = find_dimension_id_from_compare(compare)

    criteria = WmgQueryCriteria(**request["filter"])

    with ServerTiming.time("query and build response"):
        snapshot: WmgSnapshot = load_snapshot()
        q = WmgQuery(snapshot)
        default = snapshot.expression_summary_default_cube is not None and compare is None
        for dim in criteria.dict():
            if len(criteria.dict()[dim]) > 0 and depluralize(dim) in expression_summary_non_indexed_dims:
                default = False
                break

        expression_summary = q.expression_summary_default(criteria) if default else q.expression_summary(criteria)

        cell_counts = q.cell_counts(criteria)
        if expression_summary.shape[0] > 0 or cell_counts.shape[0] > 0:
            group_by_terms = ["tissue_ontology_term_id", "cell_type_ontology_term_id", compare] if compare else None

            dot_plot_matrix_df, cell_counts_cell_type_agg = get_dot_plot_data(
                expression_summary, cell_counts, group_by_terms
            )
            if is_rollup:
                dot_plot_matrix_df, cell_counts_cell_type_agg = rollup(dot_plot_matrix_df, cell_counts_cell_type_agg)

            response = jsonify(
                dict(
                    snapshot_id=snapshot.snapshot_identifier,
                    expression_summary=build_expression_summary(dot_plot_matrix_df, compare),
                    term_id_labels=dict(
                        genes=build_gene_id_label_mapping(criteria.gene_ontology_term_ids),
                        cell_types=build_ordered_cell_types_by_tissue(
                            cell_counts_cell_type_agg, snapshot.cell_type_orderings, compare
                        ),
                    ),
                )
            )
        else:  # no data, return empty json
            response = jsonify(dict(snapshot_id=snapshot.snapshot_identifier, expression_summary={}, term_id_labels={}))
    return response


def filters():
    request = connexion.request.json
    criteria = WmgFiltersQueryCriteria(**request["filter"])

    with ServerTiming.time("calculate filters and build response"):
        snapshot: WmgSnapshot = load_snapshot()
        response_filter_dims_values = build_filter_dims_values(criteria, snapshot)
        response = jsonify(
            dict(
                snapshot_id=snapshot.snapshot_identifier,
                filter_dims=response_filter_dims_values,
            )
        )
    return response


def markers():
    request = connexion.request.json
    cell_type = request["celltype"]
    tissue = request["tissue"]
    organism = request["organism"]
    n_markers = request["n_markers"]
    test = request["test"]
    snapshot: WmgSnapshot = load_snapshot()

    criteria = MarkerGeneQueryCriteria(
        tissue_ontology_term_id=tissue,
        organism_ontology_term_id=organism,
        cell_type_ontology_term_id=cell_type,
    )
    q = WmgQuery(snapshot)
    df = q.marker_genes(criteria)
    marker_genes = retrieve_top_n_markers(df, test, n_markers)
    return jsonify(
        dict(
            snapshot_id=snapshot.snapshot_identifier,
            marker_genes=marker_genes,
        )
    )


def differentialExpression():
    request = connexion.request.json

    queryGroup1Filters = request["queryGroup1Filters"]
    queryGroup2Filters = request["queryGroup2Filters"]

    criteria1 = FmgQueryCriteria(**queryGroup1Filters)
    criteria2 = FmgQueryCriteria(**queryGroup2Filters)

    snapshot: WmgSnapshot = load_snapshot()

    q = WmgQuery(snapshot)

    with ServerTiming.time("run differential expression"):
        results1, results2, success = run_differential_expression(q, criteria1, criteria2)

    return jsonify(
        dict(
            snapshot_id=snapshot.snapshot_identifier,
            differentialExpressionResults1=results1,
            differentialExpressionResults2=results2,
            tooManyCells=not success,
        )
    )


def fetch_datasets_metadata(snapshot: WmgSnapshot, dataset_ids: Iterable[str]) -> List[Dict]:
    return [
        snapshot.dataset_dict.get(dataset_id, dict(id=dataset_id, label="", collection_id="", collection_label=""))
        for dataset_id in dataset_ids
    ]


def find_dimension_id_from_compare(compare: str) -> str:
    if compare == "sex":
        return "sex_ontology_term_id"
    elif compare == "self_reported_ethnicity":
        return "self_reported_ethnicity_ontology_term_id"
    elif compare == "disease":
        return "disease_ontology_term_id"
    else:
        return None


def is_criteria_empty(criteria: WmgFiltersQueryCriteria) -> bool:
    criteria = criteria.dict()
    for key in criteria:
        if key != "organism_ontology_term_id":
            if isinstance(criteria[key], list):
                if len(criteria[key]) > 0:
                    return False
            else:
                if criteria[key] != "":
                    return False
    return True


def build_filter_dims_values(criteria: WmgFiltersQueryCriteria, snapshot: WmgSnapshot) -> Dict:
    dims = {
        "dataset_id": "",
        "disease_ontology_term_id": "",
        "sex_ontology_term_id": "",
        "development_stage_ontology_term_id": "",
        "self_reported_ethnicity_ontology_term_id": "",
        "tissue_ontology_term_id": "",
        "cell_type_ontology_term_id": "",
    }
    for dim in dims:
        dims[dim] = (
            find_all_dim_option_values(snapshot, dim)
            if is_criteria_empty(criteria)
            else find_dim_option_values(criteria, snapshot, dim)
        )

    response_filter_dims_values = dict(
        datasets=fetch_datasets_metadata(snapshot, dims["dataset_id"]),
        disease_terms=build_ontology_term_id_label_mapping(dims["disease_ontology_term_id"]),
        sex_terms=build_ontology_term_id_label_mapping(dims["sex_ontology_term_id"]),
        development_stage_terms=build_ontology_term_id_label_mapping(dims["development_stage_ontology_term_id"]),
        self_reported_ethnicity_terms=build_ontology_term_id_label_mapping(
            dims["self_reported_ethnicity_ontology_term_id"]
        ),
        tissue_terms=build_ontology_term_id_label_mapping(dims["tissue_ontology_term_id"]),
        cell_type_terms=build_ontology_term_id_label_mapping(dims["cell_type_ontology_term_id"]),
    )

    return response_filter_dims_values


def build_expression_summary(query_result: DataFrame, compare: str) -> dict:
    # Create nested dicts with gene_ontology_term_id, tissue_ontology_term_id keys, cell_type_ontology_term_id respectively
    structured_result: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(dict))
    )

    # Populate aggregated gene expressions
    query_result_agg = query_result.groupby(
        ["gene_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False
    ).agg({"nnz": "sum", "sum": "sum", "n_cells_cell_type": "sum", "n_cells_tissue": "first"})

    for i in range(query_result_agg.shape[0]):
        row = query_result_agg.iloc[i]
        structured_result[row.gene_ontology_term_id][row.tissue_ontology_term_id][row.cell_type_ontology_term_id][
            "aggregated"
        ] = dict(
            n=int(row["nnz"]),
            me=float(row["sum"] / row["nnz"]),
            pc=float(row["nnz"] / row["n_cells_cell_type"]),
            tpc=float(row["nnz"] / row["n_cells_tissue"]),
        )

    # Populate compare filter gene expressions
    if compare:
        for i in range(query_result.shape[0]):
            row = query_result.iloc[i]
            structured_result[row.gene_ontology_term_id][row.tissue_ontology_term_id][row.cell_type_ontology_term_id][
                row[compare]
            ] = dict(
                n=int(row["nnz"]),
                me=float(row["sum"] / row["nnz"]),
                pc=float(row["nnz"] / row["n_cells_cell_type"]),
                tpc=float(row["nnz"] / row["n_cells_tissue"]),
            )

    return structured_result


def agg_cell_type_counts(cell_counts: DataFrame, group_by_terms: List[str] = None) -> DataFrame:
    # Aggregate cube data by tissue, cell type
    if group_by_terms is None:
        group_by_terms = DEFAULT_GROUP_BY_TERMS
    cell_counts_cell_type_agg = cell_counts.groupby(group_by_terms, as_index=True).sum(numeric_only=True)
    cell_counts_cell_type_agg.rename(columns={"n_total_cells": "n_cells_cell_type"}, inplace=True)
    return cell_counts_cell_type_agg


def agg_tissue_counts(cell_counts: DataFrame) -> DataFrame:
    # Aggregate cube data by tissue
    cell_counts_tissue_agg = cell_counts.groupby(["tissue_ontology_term_id"], as_index=True).sum(numeric_only=True)
    cell_counts_tissue_agg.rename(columns={"n_total_cells": "n_cells_tissue"}, inplace=True)
    return cell_counts_tissue_agg


def get_dot_plot_data(
    query_result: DataFrame,
    cell_counts: DataFrame,
    group_by_terms: List[str] = None,
) -> Tuple[DataFrame, DataFrame]:
    if group_by_terms is None:
        group_by_terms = DEFAULT_GROUP_BY_TERMS
    # Get the dot plot matrix dataframe and aggregated cell counts per cell type
    cell_counts_cell_type_agg = agg_cell_type_counts(cell_counts, group_by_terms)
    cell_counts_tissue_agg = agg_tissue_counts(cell_counts)
    dot_plot_matrix_df = build_dot_plot_matrix(
        query_result, cell_counts_cell_type_agg, cell_counts_tissue_agg, group_by_terms
    )
    return dot_plot_matrix_df, cell_counts_cell_type_agg


def rollup(gene_expression_df, cell_counts_grouped_df) -> Tuple[DataFrame, DataFrame]:
    """
    Accumulates (or rolls up) cell count values and gene-expression values FOR EACH expressed gene
    up the cell type ANCESTOR paths grouped by the ontology term IDs in the multi-index of the
    input cell counts dataframe.

    Parameters
    ----------
    gene_expression_df : pandas DataFrame
        Tidy gene expression dataframe containing the dimensions across which the numeric columns will be
        aggregated.

    cell_counts_grouped_df : pandas DataFrame
        Multi-indexed cell counts dataframe containing the dimensions across which the cell count values will be
        aggregated.

    Returns
    -------
    rolled_up_gene_expression_df : pandas DataFrame
        Tidy gene expression dataframe with the same columns as the input gene expression dataframe,
        and likely greater size than the input gene expression dataframe, but with the numeric
        columns aggregated across the cell type's descendants.

    rolled_up_cell_counts_grouped_df : pandas DataFrame
        Multi-indexed cell counts dataframe with the same columns as the input cell counts dataframe,
        and likely greater size than the input cell counts dataframe, but with the cell count
        values aggregated across the cell type's descendants.
    """
    # An implementation detail note:

    # The order of operation is important: The cell counts are rolled up first to produce a the
    # `rolled_up_cell_counts_grouped_df` dataframe. Then this `rolled_up_cell_counts_grouped_df` is
    # an input into the operation to roll up gene expression values.

    # This is done for efficiency reasons where the `rolled_up_cell_counts_grouped_df` is made sparse
    # by the first operation that rolls up the cell counts. Having a sparse
    # `rolled_up_cell_counts_grouped_df` significantly improves the running time and memory footprint of
    # the second operation: rolling up the gene expression values.

    rolled_up_cell_counts_grouped_df = _rollup_cell_counts(cell_counts_grouped_df)
    rolled_up_gene_expression_df = _rollup_gene_expression(gene_expression_df, rolled_up_cell_counts_grouped_df)
    return rolled_up_gene_expression_df, rolled_up_cell_counts_grouped_df


def _add_missing_combinations_to_gene_expression_df_for_rollup(
    gene_expression_df, universal_set_cell_counts_df
) -> DataFrame:
    """
    Augments the input gene expression dataframe to include
    (gene_ontology_term_id, tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>)
    combinations for which numeric expression values should be aggregated during the rollup operation.

    Parameters
    ----------
    gene_expression_df : pandas DataFrame
        Tidy gene expression dataframe.

    universal_set_cell_counts_df : pandas DataFrame
        Multi-indexed cell counts dataframe that contains "The Universal Set Of GroupBy Ontology Terms".

    Returns
    -------
    gene_expression_with_missing_combos_df : pandas DataFrame
        Tidy gene expression dataframe with the same columns as the input gene expression dataframe,
        and likely greater size than the input gene expression dataframe that includes combinations
        for which numeric values should be aggregated during the rollup operation.
    """
    # extract group-by terms and queried genes from the input dataframes
    # if a queried gene is not present in the input dot plot dataframe, we can safely
    # ignore it as it need not be rolled up anyway.
    group_by_terms = list(universal_set_cell_counts_df.index.names)
    genes = list(set(gene_expression_df["gene_ontology_term_id"]))

    # get the names of the numeric columns
    numeric_columns = list(
        gene_expression_df.columns[[np.issubdtype(dtype, np.number) for dtype in gene_expression_df.dtypes]]
    )

    # exclude n_cells_tissue as we do not wish to roll it up
    if "n_cells_tissue" in numeric_columns:
        numeric_columns.remove("n_cells_tissue")

    # get the total number of cells per tissue to populate the n_cells_tissue in the added entries
    n_cells_tissue_dict = gene_expression_df.groupby("tissue_ontology_term_id").first()["n_cells_tissue"].to_dict()

    # get the set of available combinations of group-by terms from the aggregated cell counts
    available_combinations = set(universal_set_cell_counts_df.index.values)

    index = gene_expression_df.groupby(["gene_ontology_term_id"] + group_by_terms).first().index
    genes = index.get_level_values(0)
    combos = index.droplevel(0).values
    combinations_per_gene = to_dict(genes, combos)

    # for each gene, get the set of available combinations of group-by terms from the input expression dataframe
    entries_to_add = []
    for gene in genes:
        available_combinations_per_gene = combinations_per_gene[gene]

        # get the combinations that are missing in the input expression dataframe
        # these combinations have no data but can be rescued by the roll-up operation
        missing_combinations = available_combinations.difference(available_combinations_per_gene)
        for combo in missing_combinations:
            entry = {dim: combo[i] for i, dim in enumerate(group_by_terms)}
            entry.update({col: 0 for col in numeric_columns})
            entry["n_cells_tissue"] = n_cells_tissue_dict[entry["tissue_ontology_term_id"]]
            entry["gene_ontology_term_id"] = gene
            entries_to_add.append(entry)

    # add the missing entries to the input expression dataframe
    gene_expression_with_missing_combos_df = pd.concat((gene_expression_df, pd.DataFrame(entries_to_add)), axis=0)

    return gene_expression_with_missing_combos_df


def _build_cell_count_groups_universal_set(cell_counts_grouped_df) -> DataFrame:
    """
    Constructs a dataframe that contains all valid combination of
    (tissue, cell_type, <compare_dimension>) for which aggregation should be performed.
    We call this set of all valid combinations "The Universal Set Of GroupBy Ontology Terms".

    The combinations are encoded as a multi-index in the input cell counts dataframe.
    The "compare dimension" is optional and may not be present in the multi-index; In that case,
    the combinations considered are (tissue, cell_type) pairs.

    In this implementation, "The Universal Set Of GroupBy Ontology Terms" is constructed by
    computing the cartesian product of the (tissue, cell_type, <compare_dimension>) values
    in the input cell counts dataframe.

    Parameters
    ----------
    cell_counts_grouped_df : pandas DataFrame
        MultiIndexed cell counts dataframe containing the dimensions across which the cell count values exist.

    Returns
    -------
    universal_set_cell_counts_grouped_df : pandas DataFrame
        Multi-indexed cell counts dataframe that contains "The Universal Set Of GroupBy Ontology Terms"
        with the same columns as the input cell counts dataframe, and likely greater size than the
        input cell counts dataframe.
    """
    cartesian_product = list(
        itertools.product(
            *[cell_counts_grouped_df.index.levels[i] for i in range(len(cell_counts_grouped_df.index.names))]
        )
    )
    cartesian_product_index = pd.Index(cartesian_product)
    cartesian_product_index.set_names(cell_counts_grouped_df.index.names, inplace=True)
    universal_set_cell_counts_grouped_df = pd.DataFrame(index=cartesian_product_index)
    for c in cell_counts_grouped_df.columns:
        universal_set_cell_counts_grouped_df[c] = 0
        universal_set_cell_counts_grouped_df[c][cell_counts_grouped_df.index] = cell_counts_grouped_df[c]
    return universal_set_cell_counts_grouped_df


def _rollup_cell_counts(cell_counts_grouped_df) -> DataFrame:
    """
    Roll up cell count values across cell type descendants in the input cell counts dataframe.

    Accumulate cell count values up the cell type ontology ancestor paths for every
    (tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>) combination in the
    input cell counts dataframe. The compare dimension is an optional
    field and in the case that it is missing, then accumulation happens for every
    (tissue_ontology_term_id, cell_type_ontology_term_id) combination.

    For example:

    If (T1, C1) has cell counts in the input cell counts dataframe, and the tissue labeled T1
    has another cell labeled C2 where C2 is an ANCESTOR of C1, then the rolled up cell counts dataframe
    must include CUMULATIVE cell counts values for the combination (T1, C2) by adding in the
    cell count for (T1, C1).

    Parameters
    ----------
    cell_counts_grouped_df : pandas DataFrame
        Multi-indexed cell counts dataframe containing the dimensions across which the cell count values will be
        aggregated.

    Returns
    -------
    rolled_up_cell_counts_grouped_df : pandas DataFrame
        Multi-indexed cell counts dataframe with the same columns as the input cell counts dataframe,
        and likely greater size than the input cell counts dataframe, but with the cell count
        values aggregated across the cell type's descendants.
    """
    rolled_up_cell_counts_grouped_df = cell_counts_grouped_df

    if cell_counts_grouped_df.shape[0] > 0:
        universal_set_cell_counts_grouped_df = _build_cell_count_groups_universal_set(cell_counts_grouped_df)

        # Add index columns to the universal_set_cell_counts_grouped_df so that these columns can be
        # DIRECTLY accessed in the dataframe during the rollup operation while traversing the cell type descendants
        for col in universal_set_cell_counts_grouped_df.index.names:
            universal_set_cell_counts_grouped_df[col] = universal_set_cell_counts_grouped_df.index.get_level_values(col)

        # rollup cell counts across cell type descendants
        rolled_up_cell_counts_grouped_df = rollup_across_cell_type_descendants(universal_set_cell_counts_grouped_df)

        rolled_up_cell_counts_grouped_df = rolled_up_cell_counts_grouped_df[
            rolled_up_cell_counts_grouped_df["n_cells_cell_type"] > 0
        ]
        # Remove columns that were added to the cell counts dataframe for the purpose of rollup.
        # This is make it congruent with the structure of the input cell counts dataframe
        rolled_up_cell_counts_grouped_df.drop(columns=rolled_up_cell_counts_grouped_df.index.names, inplace=True)

    return rolled_up_cell_counts_grouped_df


def _rollup_gene_expression(gene_expression_df, universal_set_cell_counts_df) -> DataFrame:
    """
    Roll up numeric values across cell type descendants in the input gene expression dataframe.

    Accumulate gene expression values up the cell type ontology ancestor paths for every
    (gene_ontology_term_id, tissue_ontology_term_id, cell_type_ontology_term_id, <compare_dimension>)
    combination in the input gene expression dataframe. The compare dimension is an optional
    field and in the case that it is missing, then accumulation happens for every
    (gene_ontology_term_id, tissue_ontology_term_id, cell_type_ontology_term_id) combination.

    For example:

    If (G1, T1, C1) has gene expression values in the input expression dataframe,
    and the tissue labeled T1 has another cell labeled C2 where C2 is an ANCESTOR of C1,
    then the rolled up gene expression dataframe must include CUMULATIVE gene expression values for the
    combination (G1, T1, C2) by adding in the gene expression for (G1, T1, C1).

    Parameters
    ----------
    gene_expression_df : pandas DataFrame
        Tidy gene expression dataframe containing the dimensions across which the numeric columns will be
        aggregated.

    universal_set_cell_counts_df : pandas DataFrame
        Multi-indexed cell counts dataframe that contains "The Universal Set Of GroupBy Ontology Terms".

    Returns
    -------
    rolled_up_gene_expression_df : pandas DataFrame
        Tidy gene expression dataframe with the same columns as the input gene expression dataframe,
        and likely greater size than the input gene expression dataframe, but with the numeric
        columns aggregated across the cell type's descendants.
    """
    rolled_up_gene_expression_df = gene_expression_df

    if gene_expression_df.shape[0] > 0:
        # For each gene in the query, add missing combinations (tissue, cell type, compare dimension)
        # to the expression dataframe
        gene_expression_with_missing_combos_df = _add_missing_combinations_to_gene_expression_df_for_rollup(
            gene_expression_df, universal_set_cell_counts_df
        )

        # Roll up expression dataframe
        rolled_up_gene_expression_df = rollup_across_cell_type_descendants(
            gene_expression_with_missing_combos_df, ignore_cols=["n_cells_tissue"]
        )

        # Filter out the entries that were added to the dataframe that remain zero after roll-up
        rolled_up_gene_expression_df = rolled_up_gene_expression_df[rolled_up_gene_expression_df["sum"] > 0]

    return rolled_up_gene_expression_df


def build_dot_plot_matrix(
    query_result: DataFrame,
    cell_counts_cell_type_agg: DataFrame,
    cell_counts_tissue_agg: DataFrame,
    group_by_terms: List[str] = None,
) -> DataFrame:
    if group_by_terms is None:
        group_by_terms = DEFAULT_GROUP_BY_TERMS

    # Aggregate cube data by gene, tissue, cell type
    expr_summary_agg = query_result.groupby(["gene_ontology_term_id"] + group_by_terms, as_index=False).sum(
        numeric_only=True
    )
    return expr_summary_agg.join(cell_counts_cell_type_agg, on=group_by_terms, how="left").join(
        cell_counts_tissue_agg, on=["tissue_ontology_term_id"], how="left"
    )


def build_gene_id_label_mapping(gene_ontology_term_ids: List[str]) -> List[dict]:
    return [
        {gene_ontology_term_id: gene_term_label(gene_ontology_term_id)}
        for gene_ontology_term_id in gene_ontology_term_ids
    ]


def build_ontology_term_id_label_mapping(ontology_term_ids: Iterable[str]) -> List[dict]:
    return [{ontology_term_id: ontology_term_label(ontology_term_id)} for ontology_term_id in ontology_term_ids]


# getting only cell type metadata, no genes
def build_ordered_cell_types_by_tissue(
    cell_counts_cell_type_agg: DataFrame,
    cell_type_orderings: DataFrame,
    compare: str,
) -> Dict[str, Dict[str, Dict[str, Any]]]:

    distinct_tissues_cell_types = cell_counts_cell_type_agg.reset_index().rename(
        columns={"n_cells_cell_type": "n_total_cells"}
    )

    # building order for cell types for FE to use
    cell_type_orderings["order"] = range(cell_type_orderings.shape[0])

    # make a multi index
    cell_type_orderings = cell_type_orderings.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"]).first()

    indexer = list(
        zip(
            distinct_tissues_cell_types["tissue_ontology_term_id"],
            distinct_tissues_cell_types["cell_type_ontology_term_id"],
        )
    )
    indexer_bool_filter = []
    indexer_filter = []
    for index in indexer:
        indexer_bool_filter.append(index in cell_type_orderings.index)
        if index in cell_type_orderings.index:
            indexer_filter.append(index)

    joined = distinct_tissues_cell_types[indexer_bool_filter]

    for column in cell_type_orderings:
        joined[column] = list(cell_type_orderings[column][indexer_filter])

    # Remove cell types without counts
    joined = joined[joined["n_total_cells"].notnull()]

    # Create nested dicts with tissue_ontology_term_id keys, cell_type_ontology_term_id respectively
    structured_result: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(lambda: defaultdict(dict))

    # Populate aggregated gene expressions
    joined_agg = joined.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False).agg(
        {"n_total_cells": "sum", "depth": "first", "order": "first"}
    )

    agg = cell_counts_cell_type_agg.groupby(["tissue_ontology_term_id", "cell_type_ontology_term_id"]).sum().T

    for i in range(joined_agg.shape[0]):
        row = joined_agg.iloc[i]
        structured_result[row.tissue_ontology_term_id][row.cell_type_ontology_term_id]["aggregated"] = {
            "cell_type_ontology_term_id": row.cell_type_ontology_term_id,
            "name": ontology_term_label(row.cell_type_ontology_term_id),
            "total_count": int(agg[row.tissue_ontology_term_id][row.cell_type_ontology_term_id]["n_cells_cell_type"]),
            "order": int(row.order),
        }

    # Populate compare filter gene expressions
    cell_counts_cell_type_agg_T = cell_counts_cell_type_agg.T
    if compare:
        for i in range(joined.shape[0]):
            row = joined.iloc[i]
            id_to_label = build_ontology_term_id_label_mapping([row[compare]])[0]
            name = id_to_label.pop(row[compare])
            structured_result[row.tissue_ontology_term_id][row.cell_type_ontology_term_id][row[compare]] = {
                "cell_type_ontology_term_id": row.cell_type_ontology_term_id,
                "name": name if name else row[compare],
                "total_count": int(
                    cell_counts_cell_type_agg_T[row.tissue_ontology_term_id][row.cell_type_ontology_term_id][
                        row[compare]
                    ]["n_cells_cell_type"]
                ),
                "order": int(row.order),
            }

    return structured_result


def should_use_default_cube(criteria):
    default = True
    for dim in criteria.dict():
        if len(criteria.dict()[dim]) > 0 and depluralize(dim) in expression_summary_non_indexed_dims:
            default = False
            break
    return default


def de_get_expression_summary(q, cell_counts, criteria, threshold=2000):
    # all_genes = [
    #     list(i.keys())[0]
    #     for i in q._snapshot.primary_filter_dimensions["gene_terms"][criteria.organism_ontology_term_id]
    # ]

    if should_use_default_cube(criteria):
        es_agg = q.expression_summary_default(criteria).groupby("gene_ontology_term_id").sum(numeric_only=True)
    elif cell_counts.shape[0] > threshold:
        return None
    else:
        es_agg = q.expression_summary(criteria).groupby("gene_ontology_term_id").sum(numeric_only=True)

    return es_agg


def run_differential_expression(q, criteria1, criteria2, pval_thr=1e-5, threshold=2000):
    cell_counts1 = q.cell_counts(criteria1)
    cell_counts2 = q.cell_counts(criteria2)

    n_cells1 = cell_counts1["n_total_cells"].sum()
    n_cells2 = cell_counts2["n_total_cells"].sum()

    es_agg1 = de_get_expression_summary(q, cell_counts1, criteria1, threshold=threshold)
    es_agg2 = de_get_expression_summary(q, cell_counts2, criteria2, threshold=threshold)
    if es_agg1 is None or es_agg2 is None:
        return [], [], False

    genes = list(set(list(es_agg1.index) + list(es_agg2.index)))

    genes_indexer = pd.Series(index=genes, data=np.arange(len(genes)))

    sums1 = np.zeros(len(genes))
    sqsums1 = np.zeros(len(genes))

    sums2 = np.zeros(len(genes))
    sqsums2 = np.zeros(len(genes))

    sums1[genes_indexer[es_agg1.index]] = es_agg1["sum"].values
    sqsums1[genes_indexer[es_agg1.index]] = es_agg1["sqsum"].values

    sums2[genes_indexer[es_agg2.index]] = es_agg2["sum"].values
    sqsums2[genes_indexer[es_agg2.index]] = es_agg2["sqsum"].values

    pvals, effects = _run_ttest(sums1, sqsums1, n_cells1, sums2, sqsums2, n_cells2)
    de_genes = np.array(genes)[np.argsort(-effects)]
    p = pvals[np.argsort(-effects)]
    effects = effects[np.argsort(-effects)]

    statistics1 = []
    for i in range(len(p)):
        pi = p[i]
        ei = abs(effects[i])
        if ei is not np.nan and pi is not np.nan and pi < pval_thr:
            statistics1.append({"gene_ontology_term_id": de_genes[i], "p_value": pi, "effect_size": ei})
            if len(statistics1) >= 250:
                break

    pvals, effects = _run_ttest(sums2, sqsums2, n_cells2, sums1, sqsums1, n_cells1)
    de_genes = np.array(genes)[np.argsort(-effects)]
    p = pvals[np.argsort(-effects)]
    effects = effects[np.argsort(-effects)]

    statistics2 = []
    for i in range(len(p)):
        pi = p[i]
        ei = abs(effects[i])
        if ei is not np.nan and pi is not np.nan and pi < pval_thr:
            statistics2.append({"gene_ontology_term_id": de_genes[i], "p_value": pi, "effect_size": ei})
            if len(statistics2) >= 250:
                break
    return statistics1, statistics2, True
