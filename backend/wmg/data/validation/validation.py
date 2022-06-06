import logging
import os
import pathlib

import anndata
from pathlib import Path

import tiledb

from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_CUBE_NAME, CELL_COUNTS_CUBE_NAME
from backend.corpora.common.utils.math_utils import GB

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

env = os.getenv("ENV")
validation_dataset_uuid = "3de0ad6d-4378-4f62-b37b-ec0b75a50d94"
MIN_CUBE_SIZE_GB = 3 if env == "PROD" else 0.1
MIN_TISSUE_COUNT = 15 if env == "PROD" else 2
MIN_SPECIES_COUNT = 2 if env == "PROD" else 1
MIN_DATASET_COUNT = 50 if env == "PROD" else 5
MIN_MALAT1_GENE_EXPRESSION_CELL_COUNT_PERCENT = 90
MIN_ACTB_GENE_EXPRESSION_CELL_COUNT_PERCENT = 70

validation_gene_ontologies = {
    "MALAT1": "ENSG00000251562",
    "CCL5": "ENSG00000271503",
    "ACTB": "ENSG00000075624",
    "XIST": "ENSG00000229807",
    "FCN1": "ENSG00000085265",
    "TUBB4B": "ENSG00000188229",
    "CD68": "ENSG00000129226",
    "AQP5": "ENSG00000161798",
}

validation_species_ontologies = {"human": "NCBITaxon:9606", "mouse": "NCBITaxon:10090"}

validation_tissues_with_many_cell_types = {
    "lung": "UBERON:0002048",
    "blood": "UBERON:0000178",
    "lymph node": "UBERON:0000029",
    "eye": "UBERON:0000970",
    "renal medulla": "UBERON:0000362",
    "nasal cavity": "UBERON:0001707",
}
validation_cell_types = {
    "intermediate monocytes": "CL:0002393",
    "classical monocytes": "CL:0000860",
    "Non-classical monocyte": "CL:0000875",
    "ciliated cells": "CL:0000064",
    "lung ciliated cell": "CL:1000271",
    "macrophage": "CL:0000235",
    "alveolar macrophage": "CL:0000583",
    "goblet cells": "CL:0000160",
    "lung goblet cell": "CL:1000143",
    "respiratory goblet cell": "CL:0002370",
    "secretory cells": "CL:0000151",
    "mucus secreting cell": "CL:0000319",
    "serous cell of epithelium of bronchus": "CL:1000331",
}
validation_sex_ontologies = {"female": "PATO:0000383", "male": "PATO:0000384", "unknown": "unknown"}


def validate_corpus_load(anndata_object: anndata.AnnData, group_name: str, dataset_id: str):
    """
    Validate that the load looks sane
    """
    pass


def validate_cube(snapshot):
    """
    Many of these validation tests could probably be 'dryed up' or made more efficient
    Please dont do that

    The intention behind these tests is to build confidence in the validity of the cube.
    Keeping this code as clear and readable as possible (even if it is less efficient or
    repetitive) will give us more confidence in the validity of the data in the cube.
    These tests were written to be readable by a non engineer in order to be as sure as possible
    that we are validating the correct values in the cube
    """
    log_validation_details(snapshot)
    # check size
    validate_cube_size(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")

    # check species
    validate_cube_species(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")
    # todo check size of human v mouse

    # check datasets
    validate_dataset_counts(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")

    # todo list size of tissues?
    validate_tissues_in_cube(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")

    # check MALAT1 and ACTB
    validate_housekeeping_gene_expression_levels(
        f"{snapshot}/{CELL_COUNTS_CUBE_NAME}", f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}"
    )

    # check XIST appears in women but not men
    validate_sex_specific_marker_gene(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")
    # check human lung cells of particular types have marker genes
    validate_lung_cell_marker_genes(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")
    # check expression levels are correct for lung map dataset uuid 3de0ad6d-4378-4f62-b37b-ec0b75a50d94
    # genes ["MALAT1", "CCL5"]
    validate_expression_levels_for_particular_gene_dataset(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")


def log_validation_details(snapshot):
    logger.info(f"Starting cube validation for snapshot {snapshot}")
    logger.info(f"env is {env}")
    logger.info(f"MIN_CUBE_SIZE_GB is {MIN_CUBE_SIZE_GB}")
    logger.info(f"MIN_TISSUE_COUNT is {MIN_TISSUE_COUNT}")
    logger.info(f"MIN_SPECIES_COUNT is {MIN_SPECIES_COUNT}")
    logger.info(f"MIN_DATASET_COUNT is {MIN_DATASET_COUNT}")
    logger.info(f"MIN_MALAT1_GENE_EXPRESSION_CELL_COUNT_PERCENT is {MIN_MALAT1_GENE_EXPRESSION_CELL_COUNT_PERCENT}")
    logger.info(f"MIN_ACTB_GENE_EXPRESSION_CELL_COUNT_PERCENT is {MIN_ACTB_GENE_EXPRESSION_CELL_COUNT_PERCENT}")


def validate_cube_size(path_to_expression_cube):
    size_byte = sum(file.stat().st_size for file in Path(path_to_expression_cube).rglob("*"))
    size_gb = size_byte / GB
    assert size_gb > MIN_CUBE_SIZE_GB
    logger.info(f"Expression summary cube is {size_gb:.2f}GB")


def validate_cube_species(path_to_cube: str):
    with tiledb.open(path_to_cube, "r") as cube:
        species_list = cube.df[:].organism_ontology_term_id.drop_duplicates().to_list()  # todo compare speed to .unique
        species_count = len(species_list)
        assert species_count >= MIN_SPECIES_COUNT
        if env == "PROD":
            mandatory_species = validation_species_ontologies.values()
        else:
            mandatory_species = [validation_species_ontologies["human"]]
        for species in mandatory_species:
            assert species in species_list
        logger.info(f"{species_count} species included in cube")
        logger.info(f"Included species ids are: {[species for species in species_list]}")

        # todo check/log cell type per species


def validate_tissues_in_cube(path_to_cube):
    with tiledb.open(path_to_cube, "r") as cube:
        tissue_list = cube.df[:].tissue_ontology_term_id.drop_duplicates().to_list()
        tissue_count = len(tissue_list)
        assert tissue_count >= MIN_TISSUE_COUNT
        if env == "PROD":
            mandatory_tissues = validation_tissues_with_many_cell_types.values()
        else:
            mandatory_tissues = [validation_tissues_with_many_cell_types["blood"]]
        for tissue in mandatory_tissues:
            assert tissue in tissue_list
        logger.info(f"{tissue_count} tissues included in cube")
        logger.info(f"Included tissue ids are: {[tissues for tissues in tissue_list]}")

        # todo check/log cell type per tissue


def validate_housekeeping_gene_expression_levels(path_to_cell_count_cube, path_to_expression_summary):
    with tiledb.open(path_to_cell_count_cube, "r") as cell_count_cube:
        human_ontology_id = validation_species_ontologies["human"]
        cell_count_human = cell_count_cube.df[:, human_ontology_id:human_ontology_id].n_cells.sum()
        with tiledb.open(path_to_expression_summary) as cube:
            MALAT1_ont_id = validation_gene_ontologies["MALAT1"]
            MALAT1_human_expression_cube = cube.df[MALAT1_ont_id:MALAT1_ont_id, :, human_ontology_id:human_ontology_id]
            ACTB_ont_id = validation_gene_ontologies["ACTB"]
            ACTB_human_expression_cube = cube.df[ACTB_ont_id:ACTB_ont_id, :, human_ontology_id:human_ontology_id]
            MALAT1_cell_count = MALAT1_human_expression_cube.nnz.sum()
            ACTB_cell_count = ACTB_human_expression_cube.nnz.sum()
            # Most cells should express both genes, more cells should express MALAT1
            assert MALAT1_cell_count > ACTB_cell_count
            assert 100 * MALAT1_cell_count / cell_count_human > MIN_MALAT1_GENE_EXPRESSION_CELL_COUNT_PERCENT
            assert 100 * ACTB_cell_count / cell_count_human > MIN_ACTB_GENE_EXPRESSION_CELL_COUNT_PERCENT

            MALAT1_avg_expression = cube.df[MALAT1_ont_id:MALAT1_ont_id]["sum"].sum() / MALAT1_cell_count
            ACTB_avg_expression = cube.df[ACTB_ont_id:ACTB_ont_id]["sum"].sum() / ACTB_cell_count

            # parameterize high expression quartile value?
            # todo check correctness of ACTB expectation
            assert MALAT1_avg_expression > 5
            assert ACTB_avg_expression > 3


def validate_sex_specific_marker_gene(path_to_expression_summary):
    with tiledb.open(path_to_expression_summary) as cube:
        human_ontology_id = validation_species_ontologies["human"]
        sex_marker_gene_ontology_id = validation_gene_ontologies["XIST"]
        female_ontology_id = validation_sex_ontologies["female"]
        male_ontology_id = validation_sex_ontologies["male"]
        # slice cube by dimensions             gene_ontology                           organ (all)           species
        human_XIST_cube = cube.df[
            sex_marker_gene_ontology_id:sex_marker_gene_ontology_id, :, human_ontology_id:human_ontology_id
        ]

        female_xist_cube = human_XIST_cube.query(f"sex_ontology_term_id == '{female_ontology_id}'")
        male_xist_cube = human_XIST_cube.query(f"sex_ontology_term_id == '{male_ontology_id}'")

        # should be expressed in most female cells and no male cells
        assert female_xist_cube.nnz.sum() > male_xist_cube.nnz.sum()
        # should be expressed in females at a much higher rate
        female_avg_xist_expression = female_xist_cube["sum"].sum() / female_xist_cube["nnz"].sum()
        male_avg_xist_expression = male_xist_cube["sum"].sum() / male_xist_cube["nnz"].sum()
        # Todo -- why isnt male closer to 0?
        logger.info(f"female avg xist expression {female_avg_xist_expression}")
        logger.info(f"male avg xist expression {male_avg_xist_expression}")

        assert female_avg_xist_expression > male_avg_xist_expression


def validate_lung_cell_marker_genes(path_to_expression_summary):
    """
    Cells taken from human lungs have marker genes that are highly expressed in certain cell types
    gene: cell_type
    FCN1: monocytes
    TUBB4B: ciliated cells
    CD68: macrophages (alveioler)
    AQP5: goblet cells and secreting cells
    """
    human_ont_id = validation_species_ontologies["human"]
    lung_ont_id = validation_tissues_with_many_cell_types["lung"]

    # get avg expression value of gene for the celltype. That average should be greater than the avg for all
    # other cell types
    with tiledb.open(path_to_expression_summary) as cube:
        FCN1_ont_id = validation_gene_ontologies["FCN1"]
        FCN1_human_lung_cube = cube.df[FCN1_ont_id:FCN1_ont_id, lung_ont_id:lung_ont_id, human_ont_id:human_ont_id]
        validate_FCN1(FCN1_human_lung_cube)

        TUBB4B_ont_id = validation_gene_ontologies["TUBB4B"]
        TUBB4B_human_lung = cube.df[TUBB4B_ont_id:TUBB4B_ont_id, lung_ont_id:lung_ont_id, human_ont_id:human_ont_id]
        validate_TUBB4B(TUBB4B_human_lung)

        CD68_ont_id = validation_gene_ontologies["CD68"]
        CD68_human_lung = cube.df[CD68_ont_id:CD68_ont_id, lung_ont_id:lung_ont_id, human_ont_id:human_ont_id]
        validate_CD68(CD68_human_lung)

        AQP5_ont_id = validation_gene_ontologies["AQP5"]
        AQP5_human_lung = cube.df[AQP5_ont_id:AQP5_ont_id, lung_ont_id:lung_ont_id, human_ont_id:human_ont_id]
        validate_AQP5(AQP5_human_lung)


def validate_FCN1(FCN1_human_lung_cube):
    intermediate_monocyte_ontology_id = validation_cell_types["intermediate monocytes"]
    non_classical_monocyte_ontology_id = validation_cell_types["Non-classical monocyte"]
    classical_monocyte_ontology_id = validation_cell_types["classical monocytes"]
    monocyte_cell_type_ids = [
        intermediate_monocyte_ontology_id,
        non_classical_monocyte_ontology_id,
        classical_monocyte_ontology_id,
    ]

    FCN1_high_expression_cell_types = FCN1_human_lung_cube.query(
        f"cell_type_ontology_term_id in {monocyte_cell_type_ids}"
    )
    FCN1_non_high_expression_cell_types = FCN1_human_lung_cube.query(
        f"cell_type_ontology_term_id not in {monocyte_cell_type_ids}"
    )

    FCN1_high_expression_avg = (
        FCN1_high_expression_cell_types.sum()["sum"] / FCN1_high_expression_cell_types.sum()["nnz"]
    )
    FCN1_non_high_expression_avg = (
        FCN1_non_high_expression_cell_types.sum()["sum"] / FCN1_non_high_expression_cell_types.sum()["nnz"]
    )
    assert FCN1_high_expression_avg > FCN1_non_high_expression_avg


def validate_TUBB4B(TUBB4B_human_lung_cube):
    ciliated_cells_ontology_id = validation_cell_types["ciliated cells"]
    lung_ciliated_ontology_id = validation_cell_types["lung ciliated cell"]
    ciliated_cell_type_ids = [ciliated_cells_ontology_id, lung_ciliated_ontology_id]

    TUBB4B_high_expression_cell_types = TUBB4B_human_lung_cube.query(
        f"cell_type_ontology_term_id in {ciliated_cell_type_ids}"
    )
    TUBB4B_non_high_expression_cell_types = TUBB4B_human_lung_cube.query(
        f"cell_type_ontology_term_id not in {ciliated_cell_type_ids}"
    )

    TUBB4B_high_expression_avg = (
        TUBB4B_high_expression_cell_types.sum()["sum"] / TUBB4B_high_expression_cell_types.sum()["nnz"]
    )
    TUBB4B_non_high_expression_avg = (
        TUBB4B_non_high_expression_cell_types.sum()["sum"] / TUBB4B_non_high_expression_cell_types.sum()["nnz"]
    )
    assert TUBB4B_high_expression_avg > TUBB4B_non_high_expression_avg


def validate_CD68(CD68_human_lung_cube):
    macrophage_ontology_id = validation_cell_types["macrophage"]
    alvelolar_macrophage_ont_id = validation_cell_types["alveolar macrophage"]
    macrophage_cell_type_ids = [macrophage_ontology_id, alvelolar_macrophage_ont_id]

    CD68_high_expression_cell_types = CD68_human_lung_cube.query(
        f"cell_type_ontology_term_id in {macrophage_cell_type_ids}"
    )
    CD68_non_high_expression_cell_types = CD68_human_lung_cube.query(
        f"cell_type_ontology_term_id not in {macrophage_cell_type_ids}"
    )

    CD68_high_expression_avg = (
        CD68_high_expression_cell_types.sum()["sum"] / CD68_high_expression_cell_types.sum()["nnz"]
    )
    CD68_non_high_expression_avg = (
        CD68_non_high_expression_cell_types.sum()["sum"] / CD68_non_high_expression_cell_types.sum()["nnz"]
    )
    assert CD68_high_expression_avg > CD68_non_high_expression_avg


def validate_AQP5(AQP5_human_lung_cube):
    goblet_cell_ontology_id = validation_cell_types["goblet cells"]
    lung_goblet_ontology_id = validation_cell_types["lung goblet cell"]
    respiratory_goblet_ontology_id = validation_cell_types["respiratory goblet cell"]
    goblet_cell_type_ids = [goblet_cell_ontology_id, lung_goblet_ontology_id, respiratory_goblet_ontology_id]

    secretory_ontology_id = validation_cell_types["secretory cells"]
    mucus_secretory_ontology_id = validation_cell_types["mucus secreting cell"]
    bronchus_serous_ontology_id = validation_cell_types["serous cell of epithelium of bronchus"]
    secretory_cell_type_ids = [secretory_ontology_id, mucus_secretory_ontology_id, bronchus_serous_ontology_id]
    AQP5_high_expression_cell_type_ids = goblet_cell_type_ids + secretory_cell_type_ids
    AQP5_high_expression_cell_types = AQP5_human_lung_cube.query(
        f"cell_type_ontology_term_id in {AQP5_high_expression_cell_type_ids}"
    )
    AQP5_non_high_expression_cell_types = AQP5_human_lung_cube.query(
        f"cell_type_ontology_term_id not in {AQP5_high_expression_cell_type_ids}"
    )

    AQP5_high_expression_avg = (
        AQP5_high_expression_cell_types.sum()["sum"] / AQP5_high_expression_cell_types.sum()["nnz"]
    )
    AQP5_non_high_expression_avg = (
        AQP5_non_high_expression_cell_types.sum()["sum"] / AQP5_non_high_expression_cell_types.sum()["nnz"]
    )
    assert AQP5_high_expression_avg > AQP5_non_high_expression_avg


def validate_expression_levels_for_particular_gene_dataset(path_to_expression_summary):
    human_ont_id = validation_species_ontologies["human"]
    human_lung_int = validation_tissues_with_many_cell_types["lung"]
    MALAT1_ont_id = validation_gene_ontologies["MALAT1"]
    CCL5_ont_id = validation_gene_ontologies["CCL5"]
    with tiledb.open(path_to_expression_summary) as cube:
        MALAT1_human_lung_cube = cube.df[
            MALAT1_ont_id:MALAT1_ont_id, human_lung_int:human_lung_int, human_ont_id:human_ont_id
        ]
        CCL5_human_lung_cube = cube.df[
            CCL5_ont_id:CCL5_ont_id, human_lung_int:human_lung_int, human_ont_id:human_ont_id
        ]

        MALAT1_expression = MALAT1_human_lung_cube.query(f"dataset_id == '{validation_dataset_uuid}'")
        CCL5_expression = CCL5_human_lung_cube.query(f"dataset_id == '{validation_dataset_uuid}'")

        malat1_expression_sum_by_cell_type = MALAT1_expression.groupby("cell_type_ontology_term_id").sum()["sum"]
        ccl5_expression_sum_by_cell_type = CCL5_expression.groupby("cell_type_ontology_term_id").sum()["sum"]

        expected_values = anndata.read_h5ad(
            f"{pathlib.Path(__file__).parent.resolve()}/lung_map_3de0ad6d-4378-4f62-b37b-ec0b75a50d94.h5ad"
        )
        malat_expected = expected_values.obs.assign(MALAT1=expected_values.layers["rankit"].toarray()[:, 0])
        ccl5_expected = expected_values.obs.assign(CCL5=expected_values.layers["rankit"].toarray()[:, 1])

        expected_malat1_by_cell_type = malat_expected.groupby("cell_type_ontology_term_id").sum().MALAT1
        expected_ccl5_by_cell_type = ccl5_expected.groupby("cell_type_ontology_term_id").sum().CCL5
        # drop ccl5 cell types with expression value of zero (to match pipeline processing)
        expected_ccl5_by_cell_type = expected_ccl5_by_cell_type[expected_ccl5_by_cell_type != 0]

        # Todo actually compare once the rankit bug is fixed
        malat1_comparison = expected_malat1_by_cell_type.compare(malat1_expression_sum_by_cell_type)
        ccl5_comparison = expected_ccl5_by_cell_type.compare(ccl5_expression_sum_by_cell_type)
        print(malat1_comparison)
        print(ccl5_comparison)


def validate_dataset_counts(path_to_cube):
    # todo check # of datasets in dataset folder and number from relational db
    with tiledb.open(path_to_cube) as cube:
        datasets = cube.df[:].dataset_id.drop_duplicates()
        dataset_count = len(datasets)
        assert dataset_count > MIN_DATASET_COUNT
        logger.info(f"{dataset_count} datasets included in {path_to_cube}")
        logger.info(f"Included dataset ids are: {[dataset for dataset in datasets]}")
