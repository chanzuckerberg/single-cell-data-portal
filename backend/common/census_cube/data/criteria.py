from typing import List

from pydantic import BaseModel, Field


class BaseQueryCriteria(BaseModel):
    organism_ontology_term_id: str
    tissue_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    tissue_original_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    dataset_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    development_stage_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    disease_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    self_reported_ethnicity_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    sex_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)
    publication_citations: List[str] = Field(default=[], unique_items=True, min_items=0)
    cell_type_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=0)


class CensusCubeQueryCriteria(BaseQueryCriteria):
    gene_ontology_term_ids: List[str] = Field(default=[], unique_items=True, min_items=1)


class MarkerGeneQueryCriteria(BaseModel):
    organism_ontology_term_id: str  # required!
    tissue_ontology_term_id: str  # required!
    cell_type_ontology_term_id: str  # required!
