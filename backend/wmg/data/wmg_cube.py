from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME

import tiledb


def create_cube(tdb_group: str):
    """
    Create queryable cube and write to disk
    """
    pass


def create_cell_count_cube(tdb_group: str):
    """
    Create cell count cube and write to disk
    """
    uri = f"{tdb_group}/{CELL_COUNTS_CUBE_NAME}"
    with tiledb.open(f"{tdb_group}/obs") as obs:
        df = (
            obs.df[:]
            .groupby(
                by=[
                    "dataset_id",
                    "cell_type_ontology_term_id",
                    "tissue_ontology_term_id",
                    "assay_ontology_term_id",
                    "development_stage_ontology_term_id",
                    "disease_ontology_term_id",
                    "ethnicity_ontology_term_id",
                    "sex_ontology_term_id",
                    "organism_ontology_term_id",
                ],
                as_index=False,
            )
            .size()
        )

        tiledb.from_pandas(uri, df)
