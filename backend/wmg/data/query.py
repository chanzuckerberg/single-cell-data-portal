from typing import Dict, List, Union

import pandas as pd
from pandas import DataFrame
from pydantic import BaseModel, Field
from tiledb import Array

from backend.wmg.data.snapshot import WmgSnapshot


class WmgQueryCriteria(BaseModel):
    gene_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=1)  # required!
    organism_ontology_term_id: str  # required!
    tissue_ontology_term_ids: List[str] = Field(unique_items=True, min_items=1)  # required!
    tissue_original_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    dataset_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    # excluded per product requirements, but keeping in, commented-out, to reduce future head-scratching
    # assay_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    development_stage_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    disease_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    self_reported_ethnicity_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    sex_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    compare_dimension: str = Field(default=None, required=False)


class WmgQueryCriteriaV2(BaseModel):
    gene_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=1)  # required!
    organism_ontology_term_id: str  # required!
    tissue_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    tissue_original_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    dataset_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    # excluded per product requirements, but keeping in, commented-out, to reduce future head-scratching
    # assay_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    development_stage_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    disease_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    self_reported_ethnicity_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    sex_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    compare_dimension: str = Field(default=None, required=False)


class WmgFiltersQueryCriteria(BaseModel):
    organism_ontology_term_id: str  # required!
    tissue_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    tissue_original_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    dataset_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    development_stage_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    disease_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    self_reported_ethnicity_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    sex_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    cell_type_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)


class FmgQueryCriteria(BaseModel):
    organism_ontology_term_id: str  # required!
    tissue_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    # for now, we will only support finding marker genes for a single cell type.
    # this is to account for the fact that roll-up becomes much more complex when
    # multiple cell types are specified.
    cell_type_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0, max_items=1)
    tissue_original_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    dataset_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    # excluded per product requirements, but keeping in, commented-out, to reduce future head-scratching
    # assay_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    development_stage_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    disease_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    ethnicity_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    sex_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)


class MarkerGeneQueryCriteria(BaseModel):
    organism_ontology_term_id: str  # required!
    tissue_ontology_term_id: str  # required!
    cell_type_ontology_term_id: str  # required!


class WmgQuery:
    def __init__(self, snapshot: WmgSnapshot) -> None:
        self._snapshot = snapshot

    def expression_summary(self, criteria: WmgQueryCriteria) -> DataFrame:
        return self._query(
            cube=self._snapshot.expression_summary_cube,
            criteria=criteria,
            indexed_dims=[
                "gene_ontology_term_ids",
                "tissue_ontology_term_ids",
                "organism_ontology_term_id",
            ],
        )

    def expression_summary_default(self, criteria: WmgQueryCriteria) -> DataFrame:
        return self._query(
            cube=self._snapshot.expression_summary_default_cube,
            criteria=criteria,
            indexed_dims=["gene_ontology_term_ids", "tissue_ontology_term_ids", "organism_ontology_term_id"],
        )

    def expression_summary_fmg(self, criteria: FmgQueryCriteria) -> DataFrame:
        return self._query(
            cube=self._snapshot.expression_summary_fmg_cube,
            criteria=criteria,
            indexed_dims=[
                "tissue_ontology_term_ids",
                "organism_ontology_term_id",
                "cell_type_ontology_term_ids",
            ],
        )

    def marker_genes(self, criteria: MarkerGeneQueryCriteria) -> DataFrame:
        return self._query(
            cube=self._snapshot.marker_genes_cube,
            criteria=criteria,
            indexed_dims=[
                "tissue_ontology_term_ids",
                "organism_ontology_term_id",
                "cell_type_ontology_term_ids",
            ],
        )

    def cell_counts(self, criteria: Union[WmgQueryCriteria, FmgQueryCriteria]) -> DataFrame:
        cell_counts = self._query(
            cube=self._snapshot.cell_counts_cube,
            criteria=criteria.copy(exclude={"gene_ontology_term_ids"}),
            indexed_dims=["tissue_ontology_term_ids", "organism_ontology_term_id"],
        )
        cell_counts.rename(columns={"n_cells": "n_total_cells"}, inplace=True)  # expressed & non-expressed cells
        return cell_counts

    # TODO: refactor for readability: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues
    #  /chanzuckerberg/single-cell-data-portal/2133
    @staticmethod
    def _query(
        cube: Array, criteria: Union[WmgQueryCriteria, WmgQueryCriteriaV2, FmgQueryCriteria], indexed_dims: List[str]
    ) -> DataFrame:
        query_cond = ""
        attrs = {}
        for attr_name, vals in criteria.dict(exclude=set(indexed_dims)).items():
            attr = depluralize(attr_name)
            if query_cond and len(vals) > 0:
                query_cond += " and "
            if len(vals) == 1 and vals[0] != "":
                attrs[attr] = vals[0]
                query_cond += f"{attr} == val('{vals[0]}')"
            elif len(vals) > 1:
                attrs[attr] = vals
                query_cond += f"{attr} in {vals}"

        tiledb_dims_query = []
        for dim_name in indexed_dims:
            if criteria.dict()[dim_name]:
                tiledb_dims_query.append(criteria.dict()[dim_name])
            # If an "indexed" dimension is not included in the criteria,
            # then all values will be selected.
            else:
                tiledb_dims_query.append([])

        tiledb_dims_query = tuple(tiledb_dims_query)

        if isinstance(criteria, FmgQueryCriteria):
            query_result_df = pd.concat(
                cube.query(
                    cond=query_cond or None,
                    return_incomplete=True,
                    use_arrow=True,
                    attrs=["gene_ontology_term_id"],
                    dims=["tissue_ontology_term_id", "cell_type_ontology_term_id"],
                ).df[tiledb_dims_query]
            )
        elif criteria.compare_dimension is not None:
            query_result_df = pd.concat(
                cube.query(
                    cond=query_cond or None,
                    return_incomplete=True,
                    use_arrow=True,
                    attrs=["cell_type_ontology_term_id", criteria.compare_dimension],
                ).df[tiledb_dims_query]
            )
        else:
            query_result_df = pd.concat(
                cube.query(
                    cond=query_cond or None,
                    return_incomplete=True,
                    use_arrow=True,
                    attrs=["cell_type_ontology_term_id"],
                    dims=["tissue_ontology_term_id", "gene_ontology_term_id"],
                ).df[tiledb_dims_query]
            )
        return query_result_df

    def list_primary_filter_dimension_term_ids(self, primary_dim_name: str):
        return (
            self._snapshot.cell_counts_cube.query(attrs=[], dims=[primary_dim_name])
            .df[:]
            .groupby([primary_dim_name])
            .first()
            .index.tolist()
        )

    def list_grouped_primary_filter_dimensions_term_ids(
        self, primary_dim_name: str, group_by_dim: str
    ) -> Dict[str, List[str]]:
        return (
            self._snapshot.cell_counts_cube.query(attrs=[], dims=[primary_dim_name, group_by_dim])
            .df[:]
            .drop_duplicates()
            .groupby(group_by_dim)
            .agg(list)
            .to_dict()[primary_dim_name]
        )


def depluralize(attr_name):
    return attr_name[:-1]


def retrieve_top_n_markers(query_result, test, n_markers):
    """
    Retrieve the top n markers for a given cell type and test

    Arguments
    ---------
    query_result: DataFrame
        The result of a marker genes query

    test: str
        The test used to determine the top n markers. Historically, we supported both ttest
        and binomtest. Currently, only ttest is supported.

    n_markers: int
        The number of top markers to retrieve. If n_markers is 0,
        all markers are retrieved.
    """
    if test == "binomtest":
        raise ValueError("binomtest is not supported anymore")

    attrs = [f"p_value_{test}", f"effect_size_{test}"]
    col_names = ["p_value", "effect_size"]
    markers = query_result[["gene_ontology_term_id"] + attrs].rename(columns=dict(zip(attrs, col_names)))
    markers = markers[markers["effect_size"].notna()]
    if n_markers > 0:
        markers = markers.nlargest(n_markers, "effect_size")
    else:
        markers = markers.sort_values("effect_size", ascending=False)
    records = markers[["gene_ontology_term_id"] + col_names].to_dict(orient="records")
    return records
