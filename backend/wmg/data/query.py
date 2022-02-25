from typing import List

from pandas import DataFrame
from pydantic import BaseModel, Field
from tiledb import Array


class WmgQueryCriteria(BaseModel):
    gene_term_ids: list[str] = Field(default=[], unique_items=True, min_items=0)
    organism_term_id: str
    tissue_term_ids: list[str] = Field(min_items=1, unique_items=True)
    dataset_ids: list[str] = Field(default=[], unique_items=True, min_items=0)
    assay_term_ids: list[str] = Field(default=[], unique_items=True, min_items=0)
    development_stage_term_ids: list[str] = Field(default=[], unique_items=True, min_items=0)
    disease_term_ids: list[str] = Field(default=[], unique_items=True, min_items=0)
    ethnicity_term_ids: list[str] = Field(default=[], unique_items=True, min_items=0)
    sex_term_ids: list[str] = Field(default=[], unique_items=True, min_items=0)


class WmgQuery:
    def __init__(self, cube: Array) -> None:
        super().__init__()
        self._cube = cube

    def execute(self, criteria: WmgQueryCriteria) -> DataFrame:
        # Aggregate cube data by gene, tissue, cell type
        return (
            DataFrame(
                self._cube.multi_index[criteria.gene_term_ids, criteria.tissue_term_ids, criteria.organism_term_id]
            )
            .groupby(["gene_term_id", "tissue_ontology_term_id", "cell_type_ontology_term_id"])
            .sum()
        )


def build_gene_id_label_mapping(gene_term_ids) -> List[dict]:
    return [{gene_term_id: f"{gene_term_id}_label"} for gene_term_id in gene_term_ids]


def build_cell_type_id_label_mapping(cell_type_term_ids) -> List[dict]:
    return [{cell_type_term_id: f"{cell_type_term_id}_label"} for cell_type_term_id in sorted(cell_type_term_ids)]
