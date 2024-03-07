from dataclasses import dataclass


@dataclass
class ValidExplorerCxgs:
    organism_celltype_cxgs: dict[str, list[str]]
    organism_tissue_celltype_cxgs: dict[str, dict[str, list[str]]]
