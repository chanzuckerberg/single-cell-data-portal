from typing import List, Dict
from uuid import uuid4

import connexion
from flask import jsonify

DUMMY_SNAPSHOT_UUID = uuid4().hex


def primary_filter_dimensions():
    organism_terms = [dict(oid1="olbl1"), dict(oid2="olbl2")]
    tissue_type_terms = [dict(ttid1="ttlbl1"), dict(ttid2="ttlbl2")]
    result = dict(snapshot_id=DUMMY_SNAPSHOT_UUID, organism_terms=organism_terms, tissue_type_terms=tissue_type_terms)
    return jsonify(result)


def query():
    request = connexion.request.json

    gene_term_ids = request["filter"]["gene_term_ids"]
    tissue_type_term_ids = request["filter"]["tissue_type_term_ids"]

    # TODO: Move response building logic outside of API layer.
    # This currently just builds a minimal valid response with dummy data.
    cell_type_term_ids = [f"CL{n:05d}" for n in range(3)]

    def build_viz_dot(cell_type_term_id: str, me: float, pc: float, n: int, tpc: float) -> dict:
        return dict(id=cell_type_term_id, me=me, pc=pc, n=n, tpc=tpc)

    def build_viz_dots(cell_type_term_ids_: List[str]) -> List[dict]:
        return [
            build_viz_dot(cell_type_term_id, me=0.0, pc=0.0, n=0, tpc=0.0) for cell_type_term_id in cell_type_term_ids_
        ]

    def build_tissues_types(tissue_type_term_ids_: List[str]) -> Dict[str, List]:
        return dict(
            [(tissue_type_term_id, build_viz_dots(cell_type_term_ids)) for tissue_type_term_id in tissue_type_term_ids_]
        )

    def build_genes(gene_term_ids_: List[str]) -> Dict[str, dict]:
        return dict([(gene_term_id, build_tissues_types(tissue_type_term_ids)) for gene_term_id in gene_term_ids_])

    def build_gene_id_label_mapping(gene_term_ids_) -> List[dict]:
        return [{gene_term_id: f"{gene_term_id}_label"} for gene_term_id in gene_term_ids_]

    def build_cell_type_id_label_mapping(cell_type_term_ids_) -> List[dict]:
        return [{cell_type_term_id: f"{cell_type_term_id}_label"} for cell_type_term_id in cell_type_term_ids_]

    return jsonify(
        dict(
            snapshot_id=DUMMY_SNAPSHOT_UUID,
            expression_summary=build_genes(gene_term_ids),
            term_id_labels=dict(
                genes=build_gene_id_label_mapping(gene_term_ids),
                cell_types=build_cell_type_id_label_mapping(cell_type_term_ids),
            ),
        )
    )
