from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from ddtrace import tracer
from pandas import DataFrame
from tiledb import Array

from backend.common.census_cube.data.criteria import (
    BaseQueryCriteria,
    CensusCubeQueryCriteria,
    MarkerGeneQueryCriteria,
)
from backend.common.census_cube.data.schemas.cube_schema_diffexp import cell_counts_indexed_dims
from backend.common.census_cube.data.snapshot import CensusCubeSnapshot


class CensusCubeQueryParams:
    def __init__(self, cube_query_valid_attrs, cube_query_valid_dims):
        self.cube_query_valid_attrs = cube_query_valid_attrs
        self.cube_query_valid_dims = cube_query_valid_dims

    def get_indexed_dims_to_lookup_query_criteria(self, cube: Array, pluralize: bool = True) -> list[str]:
        return [self._transform_cube_index_name(i.name, pluralize) for i in cube.schema.domain]

    def get_attrs_for_cube_query(self, cube: Array) -> list[str]:
        return [i.name for i in cube.schema if i.name in self.cube_query_valid_attrs]

    def get_dims_for_cube_query(self, cube: Array) -> list[str]:
        return [i.name for i in cube.schema.domain if i.name in self.cube_query_valid_dims]

    def _transform_cube_index_name(self, index_name: str, pluralize: bool = True) -> str:
        if index_name != "organism_ontology_term_id" and pluralize:
            return index_name + "s"

        return index_name


class CensusCubeQuery:
    def __init__(self, snapshot: CensusCubeSnapshot, cube_query_params: Optional[CensusCubeQueryParams] = None) -> None:
        self._snapshot = snapshot
        self._cube_query_params = cube_query_params

    @tracer.wrap(name="expression_summary", service="wmg-api", resource="_query", span_type="wmg-api")
    def expression_summary(self, criteria: CensusCubeQueryCriteria, compare_dimension=None) -> DataFrame:
        return self._query(
            cube=self._snapshot.expression_summary_cube,
            criteria=criteria,
            compare_dimension=compare_dimension,
        )

    @tracer.wrap(name="expression_summary_default", service="wmg-api", resource="_query", span_type="wmg-api")
    def expression_summary_default(self, criteria: CensusCubeQueryCriteria) -> DataFrame:
        return self._query(
            cube=self._snapshot.expression_summary_default_cube,
            criteria=criteria,
        )

    @tracer.wrap(name="marker_genes", service="wmg-api", resource="_query", span_type="wmg-api")
    def marker_genes(self, criteria: MarkerGeneQueryCriteria) -> DataFrame:
        return self._query(
            cube=self._snapshot.marker_genes_cube,
            criteria=criteria,
        )

    @tracer.wrap(name="cell_counts", service="wmg-api", resource="_query", span_type="wmg-api")
    def cell_counts(self, criteria: BaseQueryCriteria, compare_dimension=None) -> DataFrame:
        cell_counts = self._query(
            cube=self._snapshot.cell_counts_cube,
            criteria=criteria.copy(exclude={"gene_ontology_term_ids"}),
            compare_dimension=compare_dimension,
        )
        cell_counts.rename(columns={"n_cells": "n_total_cells"}, inplace=True)  # expressed & non-expressed cells
        return cell_counts

    def cell_counts_df(self, criteria: BaseQueryCriteria) -> DataFrame:
        df = self._snapshot.cell_counts_df
        mask = np.array([True] * len(df))
        for key, values in dict(criteria).items():
            values = values if isinstance(values, list) else [values]
            key = depluralize(key)
            if key in df.columns and values:
                mask &= df[key].isin(values)
        return df[mask].rename(columns={"n_cells": "n_total_cells"})

    def cell_counts_diffexp_df(self, criteria: BaseQueryCriteria) -> DataFrame:
        df = self._snapshot.cell_counts_diffexp_df
        mask = np.array([True] * len(df))
        for key, values in dict(criteria).items():
            values = values if isinstance(values, list) else [values]
            key = depluralize(key)
            if key in df.columns and values:
                mask &= df[key].isin(values)

        return df[mask].rename(columns={"n_cells": "n_total_cells"})

    @tracer.wrap(
        name="expression_summary_and_cell_counts_diffexp", service="de-api", resource="_query", span_type="de-api"
    )
    def expression_summary_and_cell_counts_diffexp(
        self, criteria: BaseQueryCriteria, use_simple: bool
    ) -> tuple[DataFrame, DataFrame]:
        cell_counts_diffexp_df = self.cell_counts_diffexp_df(criteria)
        cell_counts_group_id_key = "group_id_simple" if use_simple else "group_id"
        cube = (
            self._snapshot.expression_summary_diffexp_simple_cube
            if use_simple
            else self._snapshot.expression_summary_diffexp_cube
        )
        group_ids = cell_counts_diffexp_df[cell_counts_group_id_key].unique().tolist()
        return (
            pd.concat(
                cube.query(
                    return_incomplete=True,
                    use_arrow=True,
                    dims=["group_id"],
                ).df[group_ids]
            ),
            cell_counts_diffexp_df,
        )

    # TODO: refactor for readability: https://app.zenhub.com/workspaces/single-cell-5e2a191dad828d52cc78b028/issues
    #  /chanzuckerberg/single-cell-data-portal/2133
    def _query(
        self,
        cube: Array,
        criteria: Union[BaseQueryCriteria, CensusCubeQueryCriteria, MarkerGeneQueryCriteria],
        compare_dimension=None,
    ) -> DataFrame:
        indexed_dims = [dim.name for dim in cube.schema.domain]

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
        criteria_dict = criteria.dict()
        for dim_name in indexed_dims:
            if dim_name not in criteria_dict:
                dim_name = pluralize(dim_name)
            if criteria_dict.get(dim_name, []):
                tiledb_dims_query.append(criteria_dict[dim_name])
            # If an "indexed" dimension is not included in the criteria,
            # then all values will be selected.
            else:
                tiledb_dims_query.append([])

        tiledb_dims_query = tuple(tiledb_dims_query)
        numeric_attrs = [attr.name for attr in cube.schema if np.issubdtype(attr.dtype, np.number)]

        # get valid attributes from schema
        # valid means it is a required column for downstream processing

        # if self._cube_query_params is None, then all attributes are valid
        attrs = self._cube_query_params.get_attrs_for_cube_query(cube) if self._cube_query_params else None
        if attrs is not None:
            if compare_dimension is not None:
                attrs.append(compare_dimension)

            attrs += numeric_attrs

        # get valid dimensions from schema

        # if self._cube_query_params is None, then all dimensions are valid
        dims = self._cube_query_params.get_dims_for_cube_query(cube) if self._cube_query_params else None

        query_result_df = pd.concat(
            cube.query(
                cond=query_cond or None,
                return_incomplete=True,
                use_arrow=True,
                attrs=attrs,
                dims=dims,
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


def should_use_simple_group_ids(criteria: BaseQueryCriteria):
    return not any(
        depluralize(key) not in cell_counts_indexed_dims and values for key, values in dict(criteria).items()
    )


def depluralize(attr_name):
    return attr_name[:-1] if attr_name[-1] == "s" else attr_name


def pluralize(attr_name):
    return attr_name + "s" if attr_name[-1] != "s" and attr_name != "organism_ontology_term_id" else attr_name


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

    attrs = ["gene_ontology_term_id", "specificity", "marker_score"]
    markers = query_result[attrs]
    markers = markers[markers["marker_score"].notna()]
    if n_markers > 0:
        markers = markers.nlargest(n_markers, "marker_score")
    else:
        markers = markers.sort_values("marker_score", ascending=False)
    records = markers[attrs].to_dict(orient="records")
    return records
