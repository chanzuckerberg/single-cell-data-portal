from typing import List

from pandas import DataFrame
from pydantic import BaseModel, Field
from tiledb import Array

from backend.de.data.schemas.expression_summary_cube_schemas import base_expression_summary_indexed_dims
from backend.de.data.snapshot import DeSnapshot
from backend.de.data.utils import depluralize, pluralize


class DeQueryCriteria(BaseModel):
    organism_ontology_term_id: str
    tissue_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    cell_type_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    dataset_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    disease_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    self_reported_ethnicity_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    sex_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)


class DeQuery:
    def __init__(self, snapshot: DeSnapshot) -> None:
        self._snapshot = snapshot

    def expression_summary(self, criteria: DeQueryCriteria) -> DataFrame:
        cardinality_per_dimension = self._snapshot.cardinality_per_dimension
        criteria_dict = criteria.dict()
        discriminatory_power = {
            dim: len(criteria_dict[dim]) / cardinality_per_dimension[dim]
            for dim in criteria_dict
            if len(criteria_dict[dim]) > 0 and dim not in base_expression_summary_indexed_dims
        }
        use_default = len(discriminatory_power) == 0

        cube_key = "default" if use_default else min(discriminatory_power, key=discriminatory_power.get)
        cube = self._snapshot.expression_summary_cubes[cube_key]

        indexed_dims = [
            pluralize(dim.name) if dim.name != "organism_ontology_term_id" else dim.name
            for dim in list(cube.schema.domain)
        ]

        return self._query(
            cube=cube,
            criteria=criteria,
            indexed_dims=indexed_dims,
        )

    def cell_counts(self, criteria: DeQueryCriteria) -> DataFrame:
        cube = self._snapshot.cell_counts_cube
        indexed_dims = [
            pluralize(dim.name) if dim.name != "organism_ontology_term_id" else dim.name
            for dim in list(cube.schema.domain)
        ]

        cell_counts = self._query(
            cube=cube,
            criteria=criteria.copy(exclude={"gene_ontology_term_ids"}),
            indexed_dims=indexed_dims,
        )
        cell_counts.rename(columns={"n_cells": "n_total_cells"}, inplace=True)  # expressed & non-expressed cells
        return cell_counts

    def _query(self, cube: Array, criteria: DeQueryCriteria, indexed_dims: List[str]) -> DataFrame:
        tiledb_dims_query = []
        for dim_name in indexed_dims:
            if criteria.dict()[dim_name]:
                tiledb_dims_query.append(criteria.dict()[dim_name])
            # If an "indexed" dimension is not included in the criteria,
            # then all values will be selected.
            else:
                tiledb_dims_query.append([])

        tiledb_dims_query = tuple(tiledb_dims_query)

        query = self._return_query(cube=cube, criteria=criteria, indexed_dims=indexed_dims)
        query_result_df = query.df[tiledb_dims_query]
        return query_result_df

    @staticmethod
    def _return_query(cube: Array, criteria: DeQueryCriteria, indexed_dims: List[str]):
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

        attr_cond = query_cond if query_cond else None
        return cube.query(cond=attr_cond, use_arrow=True)
