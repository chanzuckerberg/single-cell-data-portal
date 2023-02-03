# Examples of ad hoc queries of the WMG Expression Summary cube using the TileDB query API.

import os

import pandas as pd

from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_CUBE_NAME, _open_cube

pd.set_option("max_columns", 10)
pd.set_option("display.width", 256)

# used by _read_s3_obj() and _open_cube()
os.environ["DEPLOYMENT_STAGE"] = "dev"

if __name__ == "__main__":
    snapshot_id = "test/1662103227"  # _read_s3obj("latest_snapshot_identifier")
    cube = _open_cube(
        f's3://cellxgene-wmg-{os.environ["DEPLOYMENT_STAGE"]}/{snapshot_id}/{EXPRESSION_SUMMARY_CUBE_NAME}/'
    )

    # "dims" and "attrs" are the logical schema, other stuff is TileDB-specific config
    # "dims" are akin to indexed columns (efficiently retrieves data from disk for queried values)
    # "attrs" are akin to non-indexed columns (cannot efficiently retrieve data from disk for queried values; performs
    # filtering in caller process)
    # print(cube.schema)

    # query on gene, tissue, organism
    # dimension order matters!

    # "raw" tiledb object query (no Pandas); returns OrderDict with dict keys as column names and dict values as column
    # arrays
    print(cube["ENSG00000182149", "UBERON:0000160", :, "NCBITaxon:9606"])

    # Pandas-based tiledb query; limit to 10 rows
    print(cube.df["ENSG00000182149", "UBERON:0000160", :, "NCBITaxon:9606"][:10])

    # Query multipe value per dimension; returns cross-product of all dimension values found
    print(
        pd.DataFrame(cube.multi_index[["ENSG00000182149", "ENSG00000182150"], "UBERON:0000160", [], "NCBITaxon:9606"])
    )

    # Query, returning only selected dims and attributes
    print(
        pd.DataFrame(
            cube.query(
                dims=["tissue_ontology_term_id", "gene_ontology_term_id"], attrs=["cell_type_ontology_term_id", "sum"]
            ).df["ENSG00000182150", "UBERON:0000160", :, "NCBITaxon:9606"]
        )
    )
