import io
import logging
import os
from collections import defaultdict
from cProfile import Profile
from pstats import Stats
from typing import Any, Dict, Iterable, List

import connexion
from ddtrace import tracer
from flask import jsonify
from pandas import DataFrame
from server_timing import Timing as ServerTiming

from backend.wmg.api.common.expression_dotplot import get_dot_plot_data
from backend.wmg.api.common.rollup import rollup
from backend.wmg.api.wmg_api_config import (
    READER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
    READER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
    WMG_API_FORCE_LOAD_SNAPSHOT_ID,
    WMG_API_SNAPSHOT_SCHEMA_VERSION,
)
from backend.wmg.data.ontology_labels import gene_term_label, ontology_term_label
from backend.wmg.data.query import (
    MarkerGeneQueryCriteria,
    WmgCubeQueryParams,
    WmgFiltersQueryCriteria,
    WmgQuery,
    WmgQueryCriteriaV2,
    retrieve_top_n_markers,
)
from backend.wmg.data.schemas.cube_schema import expression_summary_non_indexed_dims
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot
from backend.wmg.data.utils import depluralize, find_all_dim_option_values, find_dim_option_values

# TODO: add cache directives: no-cache (i.e. revalidate); impl etag
#  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell-data
#  -portal/2132

logger = logging.getLogger("wmg-v2-api")


def not_cpu_time():
    times = os.times()
    return times.elapsed - (times.system + times.user)


def cprofile_query(*, is_default, q, criteria, compare=None, measure_io=False):
    profiler = Profile(not_cpu_time) if measure_io else Profile()

    if is_default:
        expression_summary = profiler.runcall(q.expression_summary_default, criteria)
        s = io.StringIO()
        stats = Stats(profiler, stream=s)
        stats.strip_dirs()
        stats.sort_stats("cumulative")
        stats.print_callers()
        logger.info(
            f"PRATHAP!!! Profiling expression_summary_default with measure_io flag set to: {measure_io}:\n{s.getvalue()}"
        )
    else:
        expression_summary = profiler.runcall(q.expression_summary, criteria, compare_dimension=compare)
        s = io.StringIO()
        stats = Stats(profiler, stream=s)
        stats.strip_dirs()
        stats.sort_stats("cumulative")
        stats.print_callers()
        logger.info(
            f"PRATHAP!!! Profiling expression_summary with measure_io flag set to: {measure_io}:\n{s.getvalue()}"
        )

    return expression_summary


@tracer.wrap(
    name="primary_filter_dimensions", service="wmg-api", resource="primary_filter_dimensions", span_type="wmg-api"
)
def primary_filter_dimensions():
    with ServerTiming.time("load snapshot"):
        snapshot: WmgSnapshot = load_snapshot(
            snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
            explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
        )

    return jsonify(snapshot.primary_filter_dimensions)


@tracer.wrap(name="query", service="wmg-api", resource="query", span_type="wmg-api")
def query():
    request = connexion.request.json
    sanitize_api_query_dict(request["filter"])

    is_rollup = request.get("is_rollup", True)
    compare = request.get("compare", None)

    if compare:
        compare = find_dimension_id_from_compare(compare)

    criteria = WmgQueryCriteriaV2(**request["filter"])

    with ServerTiming.time("load snapshot"):
        snapshot: WmgSnapshot = load_snapshot(
            snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
            explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
        )

    with ServerTiming.time("query tiledb"):
        cube_query_params = WmgCubeQueryParams(
            cube_query_valid_attrs=READER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
            cube_query_valid_dims=READER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
        )
        q = WmgQuery(snapshot, cube_query_params)
        default = snapshot.expression_summary_default_cube is not None and compare is None
        for dim in criteria.dict():
            if len(criteria.dict()[dim]) > 0 and depluralize(dim) in expression_summary_non_indexed_dims:
                default = False
                break

        expression_summary = cprofile_query(
            is_default=default, q=q, criteria=criteria, compare=compare, measure_io=True
        )

        cell_counts = q.cell_counts(criteria, compare_dimension=compare)

        # For schema-4 we filter out comma-delimited values for `self_reported_ethnicity_ontology_term_id`
        # from being included in the grouping and rollup logic per functional requirements:
        # See: https://github.com/chanzuckerberg/single-cell/issues/596
        if (compare is not None) and compare == "self_reported_ethnicity_ontology_term_id":
            expression_summary = df_not_containing_comma_delimited_ethnicity_values(expression_summary)
            cell_counts = df_not_containing_comma_delimited_ethnicity_values(cell_counts)

    with ServerTiming.time("build response"):
        if expression_summary.shape[0] > 0 or cell_counts.shape[0] > 0:
            group_by_terms = ["tissue_ontology_term_id", "cell_type_ontology_term_id", compare] if compare else None

            gene_expression_df, cell_counts_grouped_df = get_dot_plot_data(
                expression_summary, cell_counts, group_by_terms
            )
            if is_rollup:
                rolled_gene_expression_df, rolled_cell_counts_grouped_df = rollup(
                    gene_expression_df, cell_counts_grouped_df
                )

            response = jsonify(
                dict(
                    snapshot_id=snapshot.snapshot_identifier,
                    expression_summary=build_expression_summary(gene_expression_df, rolled_gene_expression_df, compare),
                    term_id_labels=dict(
                        genes=build_gene_id_label_mapping(criteria.gene_ontology_term_ids),
                        cell_types=build_ordered_cell_types_by_tissue(
                            rolled_cell_counts_grouped_df, cell_counts_grouped_df, snapshot.cell_type_orderings, compare
                        ),
                    ),
                )
            )
        else:  # no data, return empty json
            response = jsonify(dict(snapshot_id=snapshot.snapshot_identifier, expression_summary={}, term_id_labels={}))
    return response


@tracer.wrap(name="filters", service="wmg-api", resource="filters", span_type="wmg-api")
def filters():
    request = connexion.request.json
    sanitize_api_query_dict(request["filter"])

    criteria = WmgFiltersQueryCriteria(**request["filter"])

    with ServerTiming.time("load snapshot"):
        snapshot: WmgSnapshot = load_snapshot(
            snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
            explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
        )

    with ServerTiming.time("calculate filters and build response"):
        response_filter_dims_values = build_filter_dims_values(criteria, snapshot)
        response = jsonify(
            dict(
                snapshot_id=snapshot.snapshot_identifier,
                filter_dims=response_filter_dims_values,
            )
        )
    return response


@tracer.wrap(name="markers", service="wmg-api", resource="markers", span_type="wmg-api")
def markers():
    request = connexion.request.json
    cell_type = request["celltype"]
    tissue = request["tissue"]
    organism = request["organism"]
    n_markers = request["n_markers"]
    test = request["test"]
    snapshot: WmgSnapshot = load_snapshot(
        snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION,
        explicit_snapshot_id_to_load=WMG_API_FORCE_LOAD_SNAPSHOT_ID,
    )

    criteria = MarkerGeneQueryCriteria(
        tissue_ontology_term_id=tissue,
        organism_ontology_term_id=organism,
        cell_type_ontology_term_id=cell_type,
    )

    cube_query_params = WmgCubeQueryParams(
        cube_query_valid_attrs=READER_WMG_CUBE_QUERY_VALID_ATTRIBUTES,
        cube_query_valid_dims=READER_WMG_CUBE_QUERY_VALID_DIMENSIONS,
    )

    q = WmgQuery(snapshot, cube_query_params)
    df = q.marker_genes(criteria)
    marker_genes = retrieve_top_n_markers(df, test, n_markers)
    return jsonify(
        dict(
            snapshot_id=snapshot.snapshot_identifier,
            marker_genes=marker_genes,
        )
    )


def df_not_containing_comma_delimited_ethnicity_values(input_df: DataFrame) -> DataFrame:
    """
    Return a new dataframe with only the rows that DO NOT contain comma-delimited
    values in the `self_reported_ethnicity_ontology_term_id` column.

    Parameters
    ----------
    input_df: Dataframe
        A dataframe that contains `self_reported_ethnicity_ontology_term_id` column

    Returns
    -------
    A dataframe containing only the rows that do not have a comma-delimited value
    for the `self_reported_ethnicity_ontology_term_id` column
    """
    return input_df[~input_df.self_reported_ethnicity_ontology_term_id.str.contains(",")]


def sanitize_api_query_dict(query_dict: Any):
    """
    Remove invalid values in the query dictionary encoding the query API
    request body.

    The assumption is that this function is called at the beginning of the
    API function. This usage also helps mitigate query injection attacks.

    NOTE: This is a destructive operation in that it mutates `query_dict`.

    Parameters
    ----------
    query_dict : json object
        The query dictionary to sanitize.

    Returns
    -------
    None because this function mutates the function argument
    """

    # Sanitize `self_reported_ethnicity_ontology_term_ids` by removing
    # comma-delimited values because WMG does not support filtering and grouping
    # by ethnicity terms that encode mixed ethnicities encoded as a single comma-delimited string
    # value
    if "self_reported_ethnicity_ontology_term_ids" in query_dict:
        ethnicity_term_ids = query_dict["self_reported_ethnicity_ontology_term_ids"]

        ethnicity_term_ids_to_keep = [x for x in ethnicity_term_ids if "," not in x]
        query_dict["self_reported_ethnicity_ontology_term_ids"] = ethnicity_term_ids_to_keep


def fetch_datasets_metadata(snapshot: WmgSnapshot, dataset_ids: Iterable[str]) -> List[Dict]:
    return [
        snapshot.dataset_metadata.get(dataset_id, dict(id=dataset_id, label="", collection_id="", collection_label=""))
        for dataset_id in dataset_ids
    ]


def find_dimension_id_from_compare(compare: str) -> str:
    if compare == "sex":
        return "sex_ontology_term_id"
    elif compare == "self_reported_ethnicity":
        return "self_reported_ethnicity_ontology_term_id"
    elif compare == "disease":
        return "disease_ontology_term_id"
    elif compare == "publication":
        return "publication_citation"
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


@tracer.wrap(name="build_filter_dims_values", service="wmg-api", resource="filters", span_type="wmg-api")
def build_filter_dims_values(criteria: WmgFiltersQueryCriteria, snapshot: WmgSnapshot) -> Dict:
    dims = {
        "dataset_id": "",
        "disease_ontology_term_id": "",
        "sex_ontology_term_id": "",
        "development_stage_ontology_term_id": "",
        "self_reported_ethnicity_ontology_term_id": "",
        "tissue_ontology_term_id": "",
        "cell_type_ontology_term_id": "",
        "publication_citation": "",
    }
    for dim in dims:
        dims[dim] = (
            find_all_dim_option_values(snapshot, criteria.organism_ontology_term_id, dim)
            if is_criteria_empty(criteria)
            else find_dim_option_values(criteria, snapshot, dim)
        )

    # For schema-4 we filter out comma-delimited values for `self_reported_ethnicity_ontology_term_id`
    # from the options list per functional requirements:
    # See: https://github.com/chanzuckerberg/single-cell/issues/596
    ethnicity_term_ids = dims["self_reported_ethnicity_ontology_term_id"]
    dims["self_reported_ethnicity_ontology_term_id"] = [term_id for term_id in ethnicity_term_ids if "," not in term_id]

    response_filter_dims_values = dict(
        datasets=fetch_datasets_metadata(snapshot, dims["dataset_id"]),
        disease_terms=build_ontology_term_id_label_mapping(dims["disease_ontology_term_id"]),
        sex_terms=build_ontology_term_id_label_mapping(dims["sex_ontology_term_id"]),
        development_stage_terms=build_ontology_term_id_label_mapping(dims["development_stage_ontology_term_id"]),
        self_reported_ethnicity_terms=build_ontology_term_id_label_mapping(
            dims["self_reported_ethnicity_ontology_term_id"]
        ),
        publication_citations=dims["publication_citation"],
        cell_type_terms=build_ontology_term_id_label_mapping(dims["cell_type_ontology_term_id"]),
        tissue_terms=build_ontology_term_id_label_mapping(dims["tissue_ontology_term_id"]),
    )

    return response_filter_dims_values


@tracer.wrap(name="build_expression_summary", service="wmg-api", resource="query", span_type="wmg-api")
def build_expression_summary(
    unrolled_gene_expression_df: DataFrame, rolled_gene_expression_df: DataFrame, compare: str
) -> dict:
    """
    Compute and build a data structure that contains gene expression summary statistics.

    Parameters
    ----------
    unrolled_gene_expression_df: A dataframe containing unrolled gene expression values for each gene.

    rolled_gene_expression_df: A dataframe containing rolleup gene expression values.

    compare: Optional. The compare dimension to further group the gene expression summary statistics into.

    Returns
    -------
    structured_result : A nested dictionary that contains gene expression summary statistics.
    """
    # Create nested dicts with gene_ontology_term_id, tissue_ontology_term_id keys, cell_type_ontology_term_id respectively
    structured_result: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(dict))
    )

    # Populate gene expressions stats for each (gene, tissue) combination
    #
    # For aggregations that do not contain `cell_type_ontology_term_id` as a component in the compound key
    # used to group rows, we SHOULD NOT USE the gene expression dataframe that contains rolled up values.
    # That is, we should not use `rolled_gene_expression_df` for such aggregations.
    # This is because the roll up is computed over the cell_type ontology inheritance graph,
    # and such aggregations will count the rolled up values many times.
    #
    # For group-by compound keys that do not contain the `cell_type_ontology_term_id`,
    # the unrolled gene expression dataframe containing unaggregated gene expression values should be used.
    # That is, we use `unrolled_gene_expression_df` for aggregation over the group by compound key:
    # (`gene_ontology_term_id`, `tissue_ontology_term_id`)
    gene_expr_grouped_df = unrolled_gene_expression_df.groupby(
        ["gene_ontology_term_id", "tissue_ontology_term_id"], as_index=False
    ).agg({"nnz": "sum", "sum": "sum", "n_cells_tissue": "first"})

    for i in range(gene_expr_grouped_df.shape[0]):
        row = gene_expr_grouped_df.iloc[i]
        structured_result[row.gene_ontology_term_id][row.tissue_ontology_term_id]["tissue_stats"]["aggregated"] = dict(
            n=int(row["nnz"]),
            me=(float(row["sum"] / row["nnz"]) if row["nnz"] else 0.0),
            tpc=(float(row["nnz"] / row["n_cells_tissue"]) if row["n_cells_tissue"] else 0.0),
        )

    # Populate gene expressions stats for each (gene, tissue, cell_type) combination
    #
    # For aggregations that contain `cell_type_ontology_term_id` as a component in the compound key used
    # to group rows, we can use the gene expression dataframe that contains rolled up values.
    # That is, we can use `rolled_gene_expression_df` for such aggregations.
    gene_expr_grouped_df = rolled_gene_expression_df.groupby(
        ["gene_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False
    ).agg({"nnz": "sum", "sum": "sum", "n_cells_cell_type": "sum", "n_cells_tissue": "first"})

    fill_out_structured_dict_aggregated(gene_expr_grouped_df, structured_result)

    # Populate gene expression stats for each (gene, tissue, cell_type, <compare_dimension>) combination
    if compare:
        fill_out_structured_dict_compare(rolled_gene_expression_df, structured_result, compare)

    return structured_result


def fill_out_structured_dict_aggregated(gene_expr_grouped_df, structured_result):
    n = gene_expr_grouped_df["nnz"].astype("int").values
    me = (gene_expr_grouped_df["sum"] / gene_expr_grouped_df["nnz"]).values
    pc = (gene_expr_grouped_df["nnz"] / gene_expr_grouped_df["n_cells_cell_type"]).values
    tpc = (gene_expr_grouped_df["nnz"] / gene_expr_grouped_df["n_cells_tissue"]).values
    genes = gene_expr_grouped_df["gene_ontology_term_id"].values
    tissues = gene_expr_grouped_df["tissue_ontology_term_id"].values
    cell_types = gene_expr_grouped_df["cell_type_ontology_term_id"].values

    for i in range(len(n)):
        structured_result[genes[i]][tissues[i]][cell_types[i]]["aggregated"] = dict(
            n=int(n[i]),
            me=float(me[i]),
            pc=float(pc[i]),
            tpc=float(tpc[i]),
        )


def fill_out_structured_dict_compare(rolled_gene_expression_df, structured_result, compare):
    n = rolled_gene_expression_df["nnz"].astype("int").values
    me = (rolled_gene_expression_df["sum"] / rolled_gene_expression_df["nnz"]).values
    pc = (rolled_gene_expression_df["nnz"] / rolled_gene_expression_df["n_cells_cell_type"]).values
    tpc = (rolled_gene_expression_df["nnz"] / rolled_gene_expression_df["n_cells_tissue"]).values
    genes = rolled_gene_expression_df["gene_ontology_term_id"].values
    tissues = rolled_gene_expression_df["tissue_ontology_term_id"].values
    cell_types = rolled_gene_expression_df["cell_type_ontology_term_id"].values
    compare = rolled_gene_expression_df[compare].values

    for i in range(len(n)):
        structured_result[genes[i]][tissues[i]][cell_types[i]][compare[i]] = dict(
            n=int(n[i]),
            me=float(me[i]),
            pc=float(pc[i]),
            tpc=float(tpc[i]),
        )


@tracer.wrap(name="build_gene_id_label_mapping", service="wmg-api", resource="query", span_type="wmg-api")
def build_gene_id_label_mapping(gene_ontology_term_ids: List[str]) -> List[dict]:
    return [
        {gene_ontology_term_id: gene_term_label(gene_ontology_term_id)}
        for gene_ontology_term_id in gene_ontology_term_ids
    ]


def build_ontology_term_id_label_mapping(ontology_term_ids: Iterable[str]) -> List[dict]:
    return [{ontology_term_id: ontology_term_label(ontology_term_id)} for ontology_term_id in ontology_term_ids]


def fill_out_structured_tissue_agg(tissue_agg, structured_result):
    tissues = tissue_agg["tissue_ontology_term_id"].values
    n = tissue_agg["n_cells_cell_type"].values

    for i in range(len(tissues)):
        structured_result[tissues[i]]["tissue_stats"]["aggregated"] = {
            "tissue_ontology_term_id": tissues[i],
            "name": ontology_term_label(tissues[i]),
            "total_count": int(n[i]),
            "order": -1,
        }


def fill_out_structured_cell_type_agg(cell_type_agg, structured_result, ordering):
    tissues = cell_type_agg["tissue_ontology_term_id"].values
    cell_types = cell_type_agg["cell_type_ontology_term_id"].values
    n = cell_type_agg["n_cells_cell_type"].values

    for i in range(len(tissues)):
        structured_result[tissues[i]][cell_types[i]]["aggregated"] = {
            "cell_type_ontology_term_id": cell_types[i],
            "name": ontology_term_label(cell_types[i]),
            "total_count": int(n[i]),
            "order": ordering.get((tissues[i], cell_types[i]), -1),
        }


def fill_out_structured_cell_type_compare(cell_type_agg_compare, structured_result, ordering, compare):
    tissues = cell_type_agg_compare["tissue_ontology_term_id"].values
    cell_types = cell_type_agg_compare["cell_type_ontology_term_id"].values
    n = cell_type_agg_compare["n_cells_cell_type"].values
    compare = cell_type_agg_compare[compare].values

    for i in range(len(tissues)):
        id_to_label = build_ontology_term_id_label_mapping([compare[i]])[0]
        name = id_to_label.pop(compare[i])
        structured_result[tissues[i]][cell_types[i]][compare[i]] = {
            "cell_type_ontology_term_id": cell_types[i],
            "name": name if name else compare[i],
            "total_count": int(n[i]),
            "order": int(ordering.get((tissues[i], cell_types[i]), -1)),
        }


# getting only cell type metadata, no genes
@tracer.wrap(name="build_ordered_cell_types_by_tissue", service="wmg-api", resource="query", span_type="wmg-api")
def build_ordered_cell_types_by_tissue(
    rolled_cell_counts_cell_type_agg: DataFrame,
    cell_counts_cell_type_agg: DataFrame,
    cell_type_orderings: dict,
    compare: str,
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    cell_counts_cell_type_agg = cell_counts_cell_type_agg.reset_index()
    rolled_cell_counts_cell_type_agg = rolled_cell_counts_cell_type_agg.reset_index()
    # Create nested dicts with tissue_ontology_term_id keys, cell_type_ontology_term_id respectively
    structured_result: Dict[str, Dict[str, Dict[str, Any]]] = defaultdict(lambda: defaultdict(dict))

    # Populate aggregated cell counts for each tissue
    tissue_agg = cell_counts_cell_type_agg.groupby(["tissue_ontology_term_id"], as_index=False).sum(numeric_only=True)
    fill_out_structured_tissue_agg(tissue_agg, structured_result)

    # Populate aggregated cell counts for each (tissue, cell_type) combination
    cell_type_agg = rolled_cell_counts_cell_type_agg.groupby(
        ["tissue_ontology_term_id", "cell_type_ontology_term_id"], as_index=False
    ).sum(numeric_only=True)

    fill_out_structured_cell_type_agg(cell_type_agg, structured_result, cell_type_orderings)

    # Populate aggregated cell counts for each (tissue, cell_type, <compare_dimension>) combination
    if compare:
        fill_out_structured_cell_type_compare(
            rolled_cell_counts_cell_type_agg, structured_result, cell_type_orderings, compare
        )

    return structured_result
