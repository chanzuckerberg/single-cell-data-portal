from typing import List, Dict

import tiledb
from pandas import DataFrame
from pydantic import BaseModel, Field
from tiledb import Array

from backend.wmg.data.snapshot import WmgSnapshot
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


class WmgQuery:
    def __init__(self, snapshot: WmgSnapshot) -> None:
        super().__init__()
        self._snapshot = snapshot

    def expression_summary(self, criteria: WmgQueryCriteria) -> DataFrame:
        return self._query(
            self._snapshot.expression_summary_cube,
            criteria,
            indexed_dims=["gene_ontology_term_ids", "tissue_ontology_term_ids", "organism_ontology_term_id"],
        )

    def cell_counts(self, criteria: WmgQueryCriteria) -> DataFrame:
        cell_counts = self._query(
            self._snapshot.cell_counts_cube,
            criteria.copy(exclude={"gene_ontology_term_ids"}),
            indexed_dims=["tissue_ontology_term_ids", "organism_ontology_term_id"],
        )
        cell_counts.rename(columns={"n_cells": "n_total_cells"}, inplace=True)  # expressed & non-expressed cells
        return cell_counts

    # TODO: refactor for readability: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues
    #  /chanzuckerberg/single-cell-data-portal/2133
    @staticmethod
    def _query(cube: Array, criteria: WmgQueryCriteria, indexed_dims: List[str]) -> DataFrame:

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
            for attr_name, vals in criteria.dict(exclude=set(indexed_dims)).items()
            if len(vals) == 1
        }
        single_valued_attr_conds = [f"{k} == val('{v}')" for k, v in single_valued_attrs.items()]
        query_cond_expr = " and ".join(single_valued_attr_conds)
        attr_cond = tiledb.QueryCondition(query_cond_expr) if query_cond_expr else None

        tiledb_dims_query = tuple([criteria.dict()[dim_name] or EMPTY_DIM_VALUES for dim_name in indexed_dims])
        query_result_df = cube.query(attr_cond=attr_cond).df[tiledb_dims_query]

        # Filter multi-valued attribute criteria using Pandas filtering
        multi_valued_attrs = {
            depluralize(attr_name): vals
            for attr_name, vals in criteria.dict(exclude=set(indexed_dims)).items()
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
            self._snapshot.expression_summary_cube.query(attrs=[], dims=[primary_dim_name])
            .df[:]
            .groupby([primary_dim_name])
            .first()
            .index.tolist()
        )

    def list_grouped_primary_filter_dimensions_term_ids(
        self, primary_dim_name: str, group_by_dim: str
    ) -> Dict[str, List[str]]:
        # TODO: Query the cell counts cube, for efficiency:
        #  https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues/chanzuckerberg/single-cell
        #  -data-portal/2134
        return (
            self._snapshot.expression_summary_cube.query(attrs=[], dims=[primary_dim_name, group_by_dim])
            .df[:]
            .drop_duplicates()
            .groupby(group_by_dim)
            .agg(list)
            .to_dict()[primary_dim_name]
        )


def depluralize(attr_name):
    return attr_name[:-1]
