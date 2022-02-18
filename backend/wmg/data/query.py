from typing import List, Dict

DUMMY_CELL_TYPE_TERM_IDS = [f"CL{n:05d}" for n in range(3)]


def build_viz_dot(cell_type_term_id: str, me: float, pc: float, n: int, tpc: float) -> dict:
    return dict(id=cell_type_term_id, me=me, pc=pc, n=n, tpc=tpc)


def build_viz_dots(cell_type_term_ids: List[str]) -> List[dict]:
    return [build_viz_dot(cell_type_term_id, me=0.0, pc=0.0, n=0, tpc=0.0) for cell_type_term_id in cell_type_term_ids]


def build_tissues_types(tissue_term_ids: List[str], cell_type_term_ids) -> Dict[str, List]:
    return dict(
            [(tissue_term_id, build_viz_dots(cell_type_term_ids)) for tissue_term_id in tissue_term_ids]
    )


def build_genes(gene_term_ids: List[str], tissue_term_ids: List[str], cell_type_term_ids) -> Dict[str, dict]:
    return dict([(gene_term_id, build_tissues_types(tissue_term_ids, cell_type_term_ids))
                 for gene_term_id in gene_term_ids])


def build_gene_id_label_mapping(gene_term_ids) -> List[dict]:
    return [{gene_term_id: f"{gene_term_id}_label"} for gene_term_id in gene_term_ids]


def build_cell_type_id_label_mapping(cell_type_term_ids) -> List[dict]:
    return [{cell_type_term_id: f"{cell_type_term_id}_label"} for cell_type_term_id in cell_type_term_ids]


