from typing import List, Dict, Union
import tiledb
from pandas import DataFrame
from pydantic import BaseModel, Field
from tiledb import Array


class WmgQueryCriteria(BaseModel):
    gene_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=1)
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


class FmgQueryCriteria(BaseModel):
    organism_ontology_term_id: str  # required!
    tissue_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    cell_type_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
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


def expression_summary_query(cube: Array, criteria: WmgQueryCriteria) -> DataFrame:
    return _query(
        cube=cube,
        criteria=criteria,
        indexed_dims=[
            "gene_ontology_term_ids",
            "tissue_ontology_term_ids",
            "tissue_original_ontology_term_ids",
            "organism_ontology_term_id",
        ],
    )


def expression_summary_fmg_query(cube: Array, criteria: FmgQueryCriteria) -> DataFrame:
    return _query(
        cube=cube,
        criteria=criteria,
        indexed_dims=[
            "tissue_ontology_term_ids",
            "organism_ontology_term_id",
            "cell_type_ontology_term_ids",
        ],
    )


def marker_genes_query(cube: Array, criteria: MarkerGeneQueryCriteria) -> DataFrame:
    return _query(
        cube=cube,
        criteria=criteria,
        indexed_dims=[
            "tissue_ontology_term_id",
            "organism_ontology_term_id",
            "cell_type_ontology_term_id",
        ],
    )


def cell_counts_query(cube: Array, criteria: Union[WmgQueryCriteria, FmgQueryCriteria]) -> DataFrame:
    cell_counts = _query(
        cube=cube,
        criteria=criteria.copy(exclude={"gene_ontology_term_ids"}),
        indexed_dims=["tissue_ontology_term_ids", "tissue_original_ontology_term_ids", "organism_ontology_term_id"],
    )
    cell_counts.rename(columns={"n_cells": "n_total_cells"}, inplace=True)  # expressed & non-expressed cells
    return cell_counts


# TODO: refactor for readability: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues
#  /chanzuckerberg/single-cell-data-portal/2133
def _query(cube: Array, criteria: Union[WmgQueryCriteria, FmgQueryCriteria], indexed_dims: List[str]) -> DataFrame:
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

    attr_cond = tiledb.QueryCondition(query_cond) if query_cond else None

    tiledb_dims_query = []
    for dim_name in indexed_dims:
        # Don't filter on this dimension but return all "original tissues" back
        if dim_name == "tissue_original_ontology_term_ids":
            tiledb_dims_query.append([])
        elif criteria.dict()[dim_name]:
            tiledb_dims_query.append(criteria.dict()[dim_name])
        # If an "indexed" dimension is not included in the criteria,
        # then all values will be selected.
        else:
            tiledb_dims_query.append([])

    tiledb_dims_query = tuple(tiledb_dims_query)

    query_result_df = cube.query(attr_cond=attr_cond, use_arrow=True).df[tiledb_dims_query]
    return query_result_df


def list_primary_filter_dimension_term_ids(cube: Array, primary_dim_name: str):
    return cube.query(attrs=[], dims=[primary_dim_name]).df[:].groupby([primary_dim_name]).first().index.tolist()


def list_grouped_primary_filter_dimensions_term_ids(
    cube: Array, primary_dim_name: str, group_by_dim: str
) -> Dict[str, List[str]]:
    return (
        cube.query(attrs=[], dims=[primary_dim_name, group_by_dim])
        .df[:]
        .drop_duplicates()
        .groupby(group_by_dim)
        .agg(list)
        .to_dict()[primary_dim_name]
    )


def depluralize(attr_name):
    return attr_name[:-1]
