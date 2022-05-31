import logging

import anndata
from pathlib import Path

import tiledb

from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_CUBE_NAME, CELL_COUNTS_CUBE_NAME
from backend.corpora.common.utils.math_utils import GB

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

validation_dataset_uuid = "3de0ad6-4378-4f62-b37b-ec0b75a50d94"
MIN_CUBE_SIZE_GB = 3  # todo paramaterize by env
MIN_TISSUE_COUNT = 5  # todo paramaterize by env
MIN_SPECIES_COUNT = 2  # todo paramaterize by env
MIN_DATASET_COUNT = 100  # todo paramaterize by env
MIN_HOUSKEEPING_GENE_EXPRESSION_CELL_COUNT_PERCENT = 90

validation_gene_ontologies = {
    "MALAT1": "ENSG00000251562",
    "CCL5": "ENSG00000271503",
    "ACTB": "ENSG00000075624",
    "XIST": "ENSG00000229807",
    "FCN1": "ENSG00000085265",
    "TUBB4B": "ENSG00000188229",
    "CD68": "ENSG00000129226",
    "AQP5": "ENSG00000161798"
}

validation_species_ontologies = {
    "human": "NCBITaxon:9606",
    "mouse": "NCBITaxon:10090"
}

validation_tissues_with_many_cell_types = {
    "lung": "UBERON:0002048",
    "blood": "UBERON:0000178",
    "lymph node": "UBERON:0000029",
    "eye": "UBERON:0000970",
    "renal medulla": "UBERON:0000362",
    "nasal cavity": "UBERON:0001707"
}

validation_cell_types = {
    "monocytes": "CL:0002393",
    "ciliated cells": "CL:0000064",
    "macrophage": "CL:0000235",
    "goblet cells": "CL:0000160",
    "secretory cells": "CL:0000151"
}
validation_sex_ontologies = {
    "female": "PATO:0000383",
    'male': "PATO:0000384",
    "unknown": "unknown"
}


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
    repeatative) will give us more confidence in the validity of the data in the cube.
    These tests were written to be readable by a non engineer in order to be as sure as possible
    that we are validating the correct values in the cube
    """
    # check size
    validate_cube_size(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")

    # check species
    validate_cube_species(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")
    # check size of human v mouse

    # check datasets
    validate_dataset_counts(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")

    # list size of tissues?
    validate_tissues_in_cube(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")
    # check MALAT1 and ACTB
    validate_housekeeping_gene_expression_levels(f"{snapshot}/{CELL_COUNTS_CUBE_NAME}",
                                                 f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")

    # check XIST appears in women but not men
    validate_sex_specific_marker_gene(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")
    # check human lung cells of particular types have marker genes
    validate_lung_cell_marker_genes(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")
    # FCN1: monocytes
    # TUBB4B: ciliated cells, and this gene is expressed to lower values in most other cell types
    # CD68: macrophages
    # AQP5: goblet cells and secreting cells

    # check expression levels are correct for lung map dataset uuid 3de0ad6-4378-4f62-b37b-ec0b75a50d94
    # genes ["MALAT1", "CCL5"]
    validate_expression_levels_for_particular_gene_dataset(f"{snapshot}/{EXPRESSION_SUMMARY_CUBE_NAME}")


def validate_cube_size(path_to_cube):
    size_byte = sum(
        file.stat().st_size for file in Path(path_to_cube).rglob('*'))
    size_gb = size_byte / GB
    assert (size_gb > MIN_CUBE_SIZE_GB)
    logger.info(f"Expression summary cube is {size_gb}GB")


def validate_cube_species(path_to_cube: str):
    with tiledb.open(path_to_cube, "r") as cube:
        species_list = cube.df[:].organism_ontology_term_id.drop_duplicates()  # todo compare speed to .unique
        species_count = len(species_list)
        assert (species_count >= MIN_SPECIES_COUNT)
        for mandatory_species in validation_species_ontologies.keys():
            assert (mandatory_species in species_list)
        logger.info(f"{species_count} included in cube")
        logger.info(f"Included species ids are: {[species for species in species_list]}")

        # todo check/log cell type per species


def validate_tissues_in_cube(path_to_cube):
    with tiledb.open(path_to_cube, "r") as cube:
        tissue_list = cube.df[:].tissue_ontology_term_id.drop_duplicates()
        tissue_count = len(tissue_list)
        assert (tissue_count >= MIN_TISSUE_COUNT)
        for mandatory_species in validation_species_ontologies.keys():
            assert (mandatory_species in tissue_list)
        logger.info(f"{tissue_count} included in cube")
        logger.info(f"Included tissue ids are: {[tissues for tissues in tissue_list]}")

        # todo check/log cell type per tissue


def validate_housekeeping_gene_expression_levels(path_to_cell_count_cube, path_to_expression_summary):
    with tiledb.open(path_to_cell_count_cube, "r") as cell_count_cube:
        human_ontology_id = validation_species_ontologies['human']
        cell_count_human = cell_count_cube.df[:, human_ontology_id:human_ontology_id].n_cells.sum()
        with tiledb.open(path_to_expression_summary) as cube:
            human_expression_cube = cube.df[:, :, human_ontology_id:human_ontology_id]
            MALAT1_ontology_id = validation_gene_ontologies['MALAT1']
            ACTB_ontology_id = validation_gene_ontologies['XIST']
            MALAT1_cell_count = human_expression_cube.df[MALAT1_ontology_id:MALAT1_ontology_id].nnz.sum()
            ACTB_cell_count = human_expression_cube.df[ACTB_ontology_id:ACTB_ontology_id].nnz.sum()
            # Most cells should express both genes, more cells should express MALAT1
            assert (MALAT1_cell_count > ACTB_cell_count)
            assert (100 * MALAT1_cell_count / cell_count_human > MIN_HOUSKEEPING_GENE_EXPRESSION_CELL_COUNT_PERCENT)
            assert (100 * ACTB_cell_count / cell_count_human > MIN_HOUSKEEPING_GENE_EXPRESSION_CELL_COUNT_PERCENT)

            MALAT1_avg_expression = cube.df[MALAT1_ontology_id:MALAT1_ontology_id]['sum'].sum() / MALAT1_cell_count
            ACTB_avg_expression = cube.df[ACTB_ontology_id:ACTB_ontology_id]['sum'].sum() / ACTB_cell_count

            # paramaterize high expression quartile value?
            assert (MALAT1_avg_expression > 5)
            assert (ACTB_avg_expression > 3)


def validate_sex_specific_marker_gene(path_to_expression_summary):
    with tiledb.open(path_to_expression_summary) as cube:
        human_ontology_id = validation_species_ontologies['human']
        sex_marker_gene_ontology_id = validation_gene_ontologies["XIST"]
        female_ontology_id = validation_sex_ontologies['female']
        male_ontology_id = validation_sex_ontologies['male']
        # slice cube by dimensions             gene_ontology                           organ (all)           species
        human_XIST_cube = cube.df[sex_marker_gene_ontology_id:sex_marker_gene_ontology_id, :, human_ontology_id:human_ontology_id]

        female_xist_cube = human_XIST_cube.query(f"sex_ontology_term_id == \'{female_ontology_id}\'")
        male_xist_cube = human_XIST_cube.query(f"sex_ontology_term_id == \'{male_ontology_id}\'")

        # should be expressed in most female cells and no male cells
        assert (female_xist_cube.nnz.sum() > male_xist_cube.nnz.sum())
        # should be expressed in females at a much higher rate
        # todo -- do as average per cell instead of just sum
        assert (female_xist_cube['sum'].sum() > male_xist_cube['sum'].sum() * 1000)


def validate_lung_cell_marker_genes(path_to_expression_summary):
    """
    Cells taken from human lungs have marker genes that are highly expressed in certain cell types
    gene: cell_type
    FCN1: monocytes
    TUBB4B: ciliated cells, and this gene is expressed in lower values in most other cell types
    CD68: macrophages (alveioler)
    AQP5: goblet cells and secreting cells
    """
    human_ontology_id = validation_species_ontologies['human']
    lung_ontology_id = validation_tissues_with_many_cell_types['lung']

    FCN1_ontology_id = validation_gene_ontologies['FCN1']
    TUBB4B_ontology_id = validation_gene_ontologies['TUBB4B']
    CD68_ontology_id = validation_gene_ontologies['CD68']
    AQP5_ontology_id = validation_gene_ontologies['AQP5']

    monocyte_ontology_id = validation_cell_types['monocytes'] # use multiple monocyte cell types
    ciliated_cells_ontology_id = validation_cell_types['ciliated cells']
    macrophage_ontology_id = validation_cell_types['macrophage']
    goblet_cell_ontology_id = validation_cell_types['goblet cells']
    secretory_cells = validation_cell_types['secretory cells']


    # get avg expression value of gene for the celltype. That average should be greater than the avg for all other cell types (for FCN1 would use monocytes)
    with tiledb.open(path_to_expression_summary) as cube:
        human_lung_cube = cube.df[:,lung_ontology_id:lung_ontology_id, human_ontology_id:human_ontology_id]

        # todo -- what to check for??



def validate_expression_levels_for_particular_gene_dataset(path_to_expression_summary):
    human_ontology_id = validation_species_ontologies['human']
    lung_ontology_id = validation_tissues_with_many_cell_types['lung']
    MALAT1_ontology_id = validation_gene_ontologies["MALAT1"]
    CCL5_ontology_id = validation_gene_ontologies["CCL5"]
    with tiledb.open(path_to_expression_summary) as cube:
        human_lung_cube = cube.df[:,lung_ontology_id:lung_ontology_id, human_ontology_id:human_ontology_id]
        lung_df = human_lung_cube.query(f"dataset_id == '{validation_dataset_uuid}'")
        MALAT1_expression = lung_df.query(f"gene_ontology_term_id == '{MALAT1_ontology_id}'")
        CCL5_expression = lung_df.query(f"gene_ontology_term_id == '{CCL5_ontology_id}'")

        malat1_expression_sum_by_cell_type = MALAT1_expression.groupby('cell_type_ontology_term_id').sum()['sum']
        ccl5_expression_sum_by_cell_type = CCL5_expression.groupby('cell_type_ontology_term_id').sum()['sum']

        expected_values = anndata.read_h5ad('lung_map_3de0ad6d-4378-4f62-b37b-ec0b75a50d94.h5ad')
        malat_expected = expected_values.obs.assign(MALAT1=expected_values.layers['rankit'].toarray()[:, 0])
        ccl5_expected = expected_values.obs.assign(CCL5=expected_values.layers['rankit'].toarray()[:, 1])

        expected_malat1_by_cell_type = malat_expected.groupby("cell_type_ontology_term_id").sum().MALAT1
        expected_ccl5_by_cell_type = ccl5_expected.groupby("cell_type_ontology_term_id").sum().CCL5

        ## Todo actually compare once the h5ad file is swapped
        malat1_comparison = expected_malat1_by_cell_type.compare(malat1_expression_sum_by_cell_type)
        ccl5_comparison = expected_ccl5_by_cell_type.compare(ccl5_expression_sum_by_cell_type)
        print(malat1_comparison)
        print(ccl5_comparison)




def validate_dataset_counts(path_to_cube):
    # todo check # of datasets in dataset folder and number from relational db
    with tiledb.open(path_to_cube) as cube:
        datasets = cube.df[:].dataset_id.drop_duplicates()
        dataset_count = len(datasets)
        assert (dataset_count > MIN_DATASET_COUNT)
        logger.info(f"{dataset_count} datasets included in {path_to_cube}")
        logger.info(f"Included dataset ids are: {[dataset for dataset in datasets]}")
