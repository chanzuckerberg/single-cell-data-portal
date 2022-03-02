from typing import List

import tiledb
from pandas import DataFrame
from pydantic import BaseModel, Field
from tiledb import Array

UNFILTERED_SLICE = slice(None)


class WmgQueryCriteria(BaseModel):
    gene_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    organism_ontology_term_id: str  # required!
    tissue_ontology_term_ids: List[str] = Field(unique_items=True, min_items=1)  # required!
    dataset_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    assay_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    development_stage_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    disease_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    ethnicity_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    sex_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)


indexed_dims = {"gene_ontology_term_ids", "organism_ontology_term_id", "tissue_ontology_term_ids"}


class WmgQuery:
    def __init__(self, cube: Array) -> None:
        super().__init__()
        self._cube = cube

    def execute(self, criteria: WmgQueryCriteria) -> DataFrame:
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
            criteria.gene_ontology_term_ids or UNFILTERED_SLICE,
            criteria.tissue_ontology_term_ids or UNFILTERED_SLICE,
            criteria.organism_ontology_term_id or UNFILTERED_SLICE,
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

        # Aggregate cube data by gene, tissue, cell type
        return query_result_df.groupby(
            ["gene_ontology_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"]
        ).sum()


def build_gene_id_label_mapping(gene_ontology_term_ids: List[str]) -> List[dict]:
    return [
        {gene_ontology_term_id: f"{gene_ontology_term_id}_label"} for gene_ontology_term_id in gene_ontology_term_ids
    ]


def build_cell_type_id_label_mapping(cell_type_term_ids: List[str]) -> List[dict]:
    return [{cell_type_term_id: f"{cell_type_term_id}_label"} for cell_type_term_id in sorted(cell_type_term_ids)]


def depluralize(attr_name):
    return attr_name[:-1]
