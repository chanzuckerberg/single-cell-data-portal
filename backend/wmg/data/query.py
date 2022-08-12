from typing import List, Dict

import tiledb
from pandas import DataFrame
from pydantic import BaseModel, Field
from tiledb import Array

from backend.wmg.data.snapshot import WmgSnapshot

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
            cube=self._snapshot.expression_summary_cube,
            criteria=criteria,
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
        query_cond = ""
        attrs = {}
        for attr_name, vals in criteria.dict(exclude=set(indexed_dims)).items():
            attr = depluralize(attr_name)
            if query_cond and len(vals) > 0:
                query_cond += " and "
            if len(vals) == 1:
                attrs[attr] = vals[0]
                query_cond += f"{attr} == val('{vals[0]}')"
            elif len(vals) > 1:
                attrs[attr] = vals
                query_cond += f"{attr} in {vals}"

        attr_cond = tiledb.QueryCondition(query_cond) if query_cond else None

        tiledb_dims_query = tuple([criteria.dict()[dim_name] or EMPTY_DIM_VALUES for dim_name in indexed_dims])

        # FIXME: HACK of the century. Prevent realloc() error & crash when query returns an empty result. This forces
        #  two queries when there should just one.
        if (
            len(
                cube.query(attr_cond=attr_cond, attrs=attrs, dims=["organism_ontology_term_id"]).multi_index[
                    tiledb_dims_query
                ]["organism_ontology_term_id"]
            )
            == 0
        ):
            # Return an expected empty DataFrame, but without crashing, thanks to use_arrow=False
            return cube.query(attr_cond=attr_cond, use_arrow=False).df[tiledb_dims_query]

        query_result_df = cube.query(attr_cond=attr_cond, use_arrow=True).df[tiledb_dims_query]

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
