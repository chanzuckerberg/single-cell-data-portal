# Examples of ad hoc queries of the WMG Expression Summary cube using the WmgQuery class

import os

import pandas as pd

from backend.wmg.data import query
from backend.wmg.data.snapshot import (
    EXPRESSION_SUMMARY_CUBE_NAME,
    WmgSnapshot,
)
from backend.wmg.data.snapshot import _open_cube

os.environ['DEPLOYMENT_STAGE'] = 'dev'

pd.set_option('max_columns', 10)
pd.set_option('display.width', 256)


def load_snapshot(snapshot_identifier) -> WmgSnapshot:
    snapshot_base_uri = f"s3://cellxgene-wmg-dev/{snapshot_identifier}/"
    return WmgSnapshot(
        snapshot_identifier=snapshot_identifier,
        expression_summary_cube= _open_cube(f"{snapshot_base_uri}/{EXPRESSION_SUMMARY_CUBE_NAME}"),
        cell_counts_cube=None,  # _open_cube(f"{snapshot_base_uri}/{CELL_COUNTS_CUBE_NAME}"),
        cell_type_orderings=pd.DataFrame(),
        primary_filter_dimensions=pd.DataFrame()
    )


def wmg_query() -> query.WmgQuery:
    snapshot = load_snapshot("1662103227")
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
    print(
        q._snapshot.expression_summary_cube.query(use_arrow=False, attrs=["n_cells"]).df[genes, tissues, :, organism].shape
    )
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
