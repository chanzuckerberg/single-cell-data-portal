# Examples of ad hoc queries of the WMG Expression Summary cube using the WmgQuery class

import os

import pandas as pd

from backend.wmg.api import query
from backend.wmg.data.snapshot import (
    EXPRESSION_SUMMARY_CUBE_NAME,
    EXPRESSION_SUMMARY_FMG_CUBE_NAME,
    WmgSnapshot,
)
from backend.wmg.data.snapshot import _open_cube

pd.set_option("max_columns", 10)
pd.set_option("display.width", 256)

# used by _read_s3_obj() and _open_cube()
os.environ["DEPLOYMENT_STAGE"] = "dev"


def load_snapshot(snapshot_id) -> WmgSnapshot:
    cube = _open_cube(
        f's3://cellxgene-wmg-{os.environ["DEPLOYMENT_STAGE"]}/{snapshot_id}/{EXPRESSION_SUMMARY_CUBE_NAME}/'
    )
    cube_fmg = _open_cube(
        f's3://cellxgene-wmg-{os.environ["DEPLOYMENT_STAGE"]}/{snapshot_id}/{EXPRESSION_SUMMARY_FMG_CUBE_NAME}/'
    )
    return WmgSnapshot(
        snapshot_identifier=snapshot_id,
        expression_summary_cube=cube,
        cell_counts_cube=None,  # _open_cube(f"{snapshot_base_uri}/{CELL_COUNTS_CUBE_NAME}"),
        cell_type_orderings=pd.DataFrame(),
        primary_filter_dimensions=pd.DataFrame(),
        expression_summary_fmg_cube=cube_fmg,
        dataset_to_gene_ids={},
        marker_genes_cube=None,
    )


def wmg_query() -> query.WmgQuery:
    snapshot_id = 1662103227  # _read_s3obj("latest_snapshot_identifier")
    snapshot = load_snapshot(snapshot_id)
    return query.WmgQuery(snapshot)


if __name__ == "__main__":

    q = wmg_query()
    genes = ["ENSG00000182149"]
    tissues = ["UBERON:0000160"]
    # tissues = ["UBERON:0002368"]
    organism = "NCBITaxon:9606"
    # print(q._snapshot.expression_summary_cube.multi_index[genes, :, :, organism])
    #
    # print(q._snapshot.expression_summary_cube[genes[0], tissues[0], :, organism])
    # print(q._snapshot.expression_summary_cube.multi_index[genes, tissues, [], organism])
    print(q._snapshot.expression_summary_cube.query().multi_index[genes, tissues, [], organism])
    print(pd.DataFrame(q._snapshot.expression_summary_cube.multi_index[genes, tissues, [], organism]))
    print(q._snapshot.expression_summary_cube.query(attrs=["nnz"]).df[genes, tissues, :, organism].shape)
    # print(q._snapshot.expression_summary_cube.df[genes, tissues, :, organism])

    # res = q.expression_summary(
    #     query.WmgQueryCriteria(
    #         gene_ontology_term_ids=genes,
    #         organism_ontology_term_id=organism,
    #         tissue_ontology_term_ids=tissues,
    #     )
    # )
    # print(res)

    # res = q.cell_counts(
    #     query.WmgQueryCriteria(
    #         organism_ontology_term_id=organism,
    #         tissue_ontology_term_ids=tissues,
    #     )
    # )
    # print(res)
