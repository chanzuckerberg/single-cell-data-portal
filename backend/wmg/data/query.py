from typing import List, Dict

import tiledb
from pandas import DataFrame
from pydantic import BaseModel, Field
from tiledb import Array

from backend.wmg.data.ontology_labels import gene_term_label, ontology_term_label

EMPTY_DIM_VALUES = ""


class WmgQueryCriteria(BaseModel):
    gene_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=1)
    organism_ontology_term_id: str  # required!
    tissue_ontology_term_ids: List[str] = Field(unique_items=True, min_items=1)  # required!
    dataset_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    # excluded per product requirements, but keeping in, commented-out, to reduce future head-scratching
    # assay_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    development_stage_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    disease_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    ethnicity_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    sex_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)


indexed_dims = {"gene_ontology_term_ids", "organism_ontology_term_id", "tissue_ontology_term_ids"}


class WmgQuery:
    def __init__(self, cube: Array) -> None:
        super().__init__()
        self._cube = cube

    # TODO: refactor for readability: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues
    #  /chanzuckerberg/single-cell-data-portal/2133
    def expression_summary(self, criteria: WmgQueryCriteria) -> DataFrame:
        # As TileDB API does not yet support logical OR'ing of attribute values in query conditions, the best we can do
        # for a single TileDB query is to have TileDB perform the filtering for only the attributes that have
        # a single criterion value specified by the user. We can then perform the filtering client-side for the
        # multi-valued attributes, albeit at the cost of retrieving extra data from the TileDB engine. Note, however,
        # that this extra data is always being retrieved by TileDB from the filesystem (and so across the network), so
        # the performance impact should not be as bad as it may at first appear.

        # Filter single-valued attribute criterion using TileDB filtering, using QueryCondition
        # TODO: Benchmark whether this is faster than just doing all the attribute filtering in Pandas
        single_valued_attrs = {
            depluralize(attr_name): vals[0]
            for attr_name, vals in criteria.dict(exclude=indexed_dims).items()
            if len(vals) == 1
        }
        single_valued_attr_conds = [f"{k} == val('{v}')" for k, v in single_valued_attrs.items()]
        query_cond_expr = " and ".join(single_valued_attr_conds)
        attr_cond = tiledb.QueryCondition(query_cond_expr) if query_cond_expr else None

        tiledb_dims_query = (
            # use of EMPTY_DIM_VALUES ensures an empty result if any of the primary dimensions are unspecified (
            # avoids introducing another code path for forming an empty result return value)
            criteria.gene_ontology_term_ids or EMPTY_DIM_VALUES,
            criteria.tissue_ontology_term_ids or EMPTY_DIM_VALUES,
            criteria.organism_ontology_term_id or EMPTY_DIM_VALUES,
        )
        query_result_df = self._cube.query(attr_cond=attr_cond).df[tiledb_dims_query]

        # Filter multi-valued attribute criteria using Pandas filtering
        multi_valued_attrs = {
            depluralize(attr_name): vals
            for attr_name, vals in criteria.dict(exclude=indexed_dims).items()
            if len(vals) > 1
        }
        if multi_valued_attrs:
            # Note that this is a Pandas query expression, which happens to be very similar to TileDB Python
            # QueryCondition expression. But don't confuse the two! Pandas query expressions _do_ support logical OR'ing
            # in filtering expressions
            multi_valued_attr_conds = [f"{k} in {v}" for k, v in multi_valued_attrs.items()]
            query_result_df = query_result_df.query(" and ".join(multi_valued_attr_conds))

        return query_result_df

    def list_primary_filter_dimension_term_ids(self, primary_dim_name: str):
        # TODO: Query the cell counts cube, for efficiency:
        #  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell
        #  -data-portal/2134
        return (
            self._cube.query(attrs=[], dims=[primary_dim_name]).df[:].groupby([primary_dim_name]).first().index.tolist()
        )

    def list_grouped_primary_filter_dimensions_term_ids(
        self, primary_dim_name: str, group_by_dim: str
    ) -> Dict[str, List[str]]:
        # TODO: Query the cell counts cube, for efficiency:
        #  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell
        #  -data-portal/2134
        return (
            self._cube.query(attrs=[], dims=[primary_dim_name, group_by_dim])
            .df[:]
            .drop_duplicates()
            .groupby(group_by_dim)
            .agg(list)
            .to_dict()[primary_dim_name]
        )


def build_dot_plot_matrix(query_result: DataFrame) -> DataFrame:
    # Aggregate cube data by gene, tissue, cell type
    return query_result.groupby(
        ["gene_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"]
    ).sum()


def build_gene_id_label_mapping(gene_ontology_term_ids: List[str]) -> List[dict]:
    return [
        {gene_ontology_term_id: gene_term_label(gene_ontology_term_id)}
        for gene_ontology_term_id in gene_ontology_term_ids
    ]


def build_ontology_term_id_label_mapping(cell_type_term_ids: List[str]) -> List[dict]:
    return [
        {cell_type_term_id: ontology_term_label(cell_type_term_id)} for cell_type_term_id in sorted(cell_type_term_ids)
    ]


def depluralize(attr_name):
    return attr_name[:-1]
