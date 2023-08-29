import logging
import os
import pathlib
from pathlib import Path

import anndata
import tiledb

from backend.common.utils.math_utils import GB
from backend.wmg.data.snapshot import CELL_COUNTS_CUBE_NAME, EXPRESSION_SUMMARY_CUBE_NAME
from backend.wmg.data.validation import fixtures

logger: logging.Logger = logging.getLogger(__name__)


class Validation:
    def __init__(self, corpus_path):
        self.errors = []
        self.corpus_path = corpus_path
        self.expression_summary_path = f"{corpus_path}/{EXPRESSION_SUMMARY_CUBE_NAME}"
        self.cell_count_path = f"{corpus_path}/{CELL_COUNTS_CUBE_NAME}"
        self.env = os.getenv("DEPLOYMENT_STAGE")
        self.validation_dataset_id = "3de0ad6d-4378-4f62-b37b-ec0b75a50d94"
        self.MIN_CUBE_SIZE_GB = 2.3
        self.MIN_TISSUE_COUNT = 53
        self.MIN_SPECIES_COUNT = 4
        self.MIN_DATASET_COUNT = 254
        self.MIN_MALAT1_GENE_EXPRESSION_CELL_COUNT_PERCENT = 80
        self.MIN_ACTB_GENE_EXPRESSION_CELL_COUNT_PERCENT = 60
        self.MIN_MALAT1_RANKIT_EXPRESSION = 4
        self.MIN_ACTB_RANKIT_EXPRESSION = 2.75

    def validate_cube(self):
        """
        Many of these validation tests could probably be 'dryed up' or made more efficient
        Please dont do that

        The intention behind these tests is to build confidence in the validity of the cube.
        Keeping this code as clear and readable as possible (even if it is less efficient or
        repetitive) will give us more confidence in the validity of the data in the cube.
        These tests were written to be readable by a non engineer in order to be as sure as possible
        that we are validating the correct values in the cube
        """
        self.log_validation_details()

        # check size
        self.validate_cube_size()

        # check species
        self.validate_cube_species()
        # todo check size of human v mouse

        # check datasets
        self.validate_dataset_counts()

        # todo list size of tissues?
        self.validate_tissues_in_cube()

        # check tissue roll up
        self.validate_tissue_rollup_cell_count()
        self.validate_tissue_rollup_expression()

        # check MALAT1 and ACTB
        self.validate_housekeeping_gene_expression_levels()

        # check XIST appears in women but not men
        self.validate_sex_specific_marker_gene()
        # check human lung cells of particular types have marker genes
        self.validate_lung_cell_marker_genes()
        # check expression levels are correct for lung map dataset uuid 3de0ad6d-4378-4f62-b37b-ec0b75a50d94
        # genes ["MALAT1", "CCL5"]
        self.validate_expression_levels_for_particular_gene_dataset()

        if len(self.errors) > 0:
            error_message = f"Cube Validation Failed with {len(self.errors)} errors"
            for error in self.errors:
                error_message += f"\n{error}"
            logger.error(error_message)
            return False
        return True

    def log_validation_details(self):
        logger.info(f"Starting cube validation for snapshot {self.corpus_path}")
        logger.info(f"env is {self.env}")
        logger.info(f"MIN_CUBE_SIZE_GB is {self.MIN_CUBE_SIZE_GB}")
        logger.info(f"MIN_TISSUE_COUNT is {self.MIN_TISSUE_COUNT}")
        logger.info(f"MIN_SPECIES_COUNT is {self.MIN_SPECIES_COUNT}")
        logger.info(f"MIN_DATASET_COUNT is {self.MIN_DATASET_COUNT}")
        logger.info(
            f"MIN_MALAT1_GENE_EXPRESSION_CELL_COUNT_PERCENT is {self.MIN_MALAT1_GENE_EXPRESSION_CELL_COUNT_PERCENT}"
        )
        logger.info(
            f"MIN_ACTB_GENE_EXPRESSION_CELL_COUNT_PERCENT is {self.MIN_ACTB_GENE_EXPRESSION_CELL_COUNT_PERCENT}"
        )
        logger.info(f"MIN_MALAT1_RANKIT_EXPRESSION is {self.MIN_MALAT1_RANKIT_EXPRESSION}")
        logger.info(f"MIN_ACTB_RANKIT_EXPRESSION is {self.MIN_ACTB_RANKIT_EXPRESSION}")

    def validate_cube_size(self):
        size_byte = sum(file.stat().st_size for file in Path(self.expression_summary_path).rglob("*"))
        size_gb = size_byte / GB
        logger.info(f"Expression summary cube is {size_gb:.2f}GB")
        if not size_gb > self.MIN_CUBE_SIZE_GB:
            self.errors.append(
                f"Expression summary cube is {size_gb:.2f}GB which is smaller than expected {self.MIN_CUBE_SIZE_GB}GB"
            )

    def validate_cube_species(self):
        with tiledb.open(self.cell_count_path, "r") as cube:
            species_list = cube.df[:].organism_ontology_term_id.drop_duplicates().to_list()
            species_count = len(species_list)
            if species_count < self.MIN_SPECIES_COUNT:
                self.errors.append(
                    f"Expression summary cube missing mandatory species. Only contains {species_count} species"
                )
            for species in fixtures.validation_species_ontologies.values():
                if species not in species_list:
                    self.errors.append(f"Cube missing species: {species}")
            logger.info(f"{species_count} species included in cube")
            logger.info(f"Included species ids are: {species_list}")

            # todo check/log cell type per species

    def validate_tissue_rollup_cell_count(self):
        """
        Validates that tissues rolled up in the axis "tissue_ontology_term_id" was correctly done from
        attribute "tissue_original_ontology_term_id".

        The way this validation works is by creating two cell count tables:
            1. Create cell count tables for Lung using "tissue_ontology_term_id"
            2. Create cell count tables for all Lung parts using "tissue_original_ontology_term_id"

        The two should be identical
        """

        with tiledb.open(self.cell_count_path, "r") as cell_count_cube:
            human_ontology_id = fixtures.validation_species_ontologies["human"]
            all_lung_tissues = fixtures.validation_all_lung_tissues
            lung_high_level_tissue = fixtures.validation_lung_high_level
            criteria_original = dict(
                tissue_original_ontology_term_id=all_lung_tissues,
                organism_ontology_term_id=human_ontology_id,
            )
            criteria_rollup = dict(
                tissue_ontology_term_id=lung_high_level_tissue,
                organism_ontology_term_id=human_ontology_id,
            )
            # the cell counts cube is small enough to fully load into memory
            cell_counts_df = cell_count_cube.df[:]
            original_tissues_cell_count = query_pandas(cell_counts_df, criteria_original).n_cells.sum()
            rollup_tissues_cell_count = query_pandas(cell_counts_df, criteria_rollup).n_cells.sum()

            if original_tissues_cell_count != rollup_tissues_cell_count:
                logger.error(
                    f"Tissue roll up error, cell counts for lung subparts ({original_tissues_cell_count}) "
                    f"is not equal to cell counts for rolled-up lung ({rollup_tissues_cell_count})"
                )

    def validate_tissue_rollup_expression(self):
        """
        Validates that tissues rolled up in the axis "tissue_ontology_term_id" was correctly done from
        attribute "tissue_original_ontology_term_id".

        The way this validation works is by creating two cell type EXPRESISION tables:
            1. Create cell type EXPRESSION tables for Lung using "tissue_ontology_term_id"
            2. Create cell type EXPRESSION tables for all Lung parts using "tissue_original_ontology_term_id"

        Expression is checked on MALAT1 gene only

        The two should be identical
        """

        with tiledb.open(self.expression_summary_path, "r") as cube:
            human_ontology_id = fixtures.validation_species_ontologies["human"]
            MALAT1_ont_id = fixtures.validation_gene_ontologies["MALAT1"]
            all_lung_tissues = fixtures.validation_all_lung_tissues
            lung_high_level_tissue = fixtures.validation_lung_high_level

            # expression_summary is large enough that we don't want to load it into memory
            # query the cube for the gene and organism and then filter the resulting dataframe
            criteria_original = dict(
                tissue_original_ontology_term_id=all_lung_tissues,
            )
            criteria_rollup = dict(
                tissue_ontology_term_id=lung_high_level_tissue,
            )

            original_tissues_expression = cube.df[MALAT1_ont_id, :, human_ontology_id]
            rollup_tissues_expression = cube.df[MALAT1_ont_id, :, human_ontology_id]

            original_tissues_expression = query_pandas(original_tissues_expression, criteria_original)
            rollup_tissues_expression = query_pandas(rollup_tissues_expression, criteria_rollup)

            original_tissues_expression = original_tissues_expression[["cell_type_ontology_term_id", "sum"]]
            rollup_tissues_expression = rollup_tissues_expression[["cell_type_ontology_term_id", "sum"]]

            original_tissues_expression = original_tissues_expression.groupby("cell_type_ontology_term_id").sum(
                numeric_only=True
            )
            rollup_tissues_expression = rollup_tissues_expression.groupby("cell_type_ontology_term_id").sum(
                numeric_only=True
            )

            if not original_tissues_expression["sum"].equals(rollup_tissues_expression["sum"]):
                logger.error(
                    f"Tissue roll up error, cell expresion for lung subparts ({original_tissues_expression}) "
                    f"is not equal to expression for rolled-up lung ({rollup_tissues_expression})"
                )

    def validate_tissues_in_cube(self):
        with tiledb.open(self.cell_count_path, "r") as cube:
            tissue_list = cube.df[:].tissue_ontology_term_id.drop_duplicates().to_list()
            tissue_count = len(tissue_list)
            if tissue_count < self.MIN_TISSUE_COUNT:
                self.errors.append(f"Only {tissue_count} tissues included in cube")
            for tissue in fixtures.validation_tissues_with_many_cell_types.values():
                if tissue not in tissue_list:
                    self.errors.append(f"{tissue} missing from tissue list")
            logger.info(f"{tissue_count} tissues included in cube")
            logger.info(f"Included tissue ids are: {tissue_list}")

            # todo check/log cell type per tissue

    def validate_housekeeping_gene_expression_levels(self):
        with tiledb.open(self.cell_count_path, "r") as cell_count_cube:
            human_ontology_id = fixtures.validation_species_ontologies["human"]
            cell_count_human = cell_count_cube.df[:, human_ontology_id].n_cells.sum()
            with tiledb.open(self.expression_summary_path) as cube:
                MALAT1_ont_id = fixtures.validation_gene_ontologies["MALAT1"]
                MALAT1_human_expression_cube = cube.df[MALAT1_ont_id, :, human_ontology_id]
                ACTB_ont_id = fixtures.validation_gene_ontologies["ACTB"]
                ACTB_human_expression_cube = cube.df[ACTB_ont_id, :, human_ontology_id]
                MALAT1_cell_count = MALAT1_human_expression_cube.nnz.sum()
                ACTB_cell_count = ACTB_human_expression_cube.nnz.sum()
                # Most cells should express both genes, more cells should express MALAT1
                if ACTB_cell_count > MALAT1_cell_count:
                    self.errors.append(f"More cells express ACTB ({ACTB_cell_count}) than MALAT1 ({MALAT1_cell_count})")
                if 100 * MALAT1_cell_count / cell_count_human < self.MIN_MALAT1_GENE_EXPRESSION_CELL_COUNT_PERCENT:
                    self.errors.append(
                        f"less than " f"{self.MIN_MALAT1_GENE_EXPRESSION_CELL_COUNT_PERCENT}% of cells express MALAT1"
                    )
                if 100 * ACTB_cell_count / cell_count_human < self.MIN_ACTB_GENE_EXPRESSION_CELL_COUNT_PERCENT:
                    self.errors.append(
                        f"less than " f"{self.MIN_ACTB_GENE_EXPRESSION_CELL_COUNT_PERCENT}% of cells express ACTB"
                    )

                MALAT1_avg_expression = cube.df[MALAT1_ont_id]["sum"].sum() / MALAT1_cell_count
                ACTB_avg_expression = cube.df[ACTB_ont_id]["sum"].sum() / ACTB_cell_count
                if MALAT1_avg_expression < self.MIN_MALAT1_RANKIT_EXPRESSION:
                    self.errors.append(f"MALAT1 avg rankit score is {MALAT1_avg_expression}")
                if ACTB_avg_expression < self.MIN_ACTB_RANKIT_EXPRESSION:
                    self.errors.append(f"ACTB avg rankit score is {ACTB_avg_expression}")

    def validate_sex_specific_marker_gene(self):
        with tiledb.open(self.expression_summary_path) as cube:
            human_ontology_id = fixtures.validation_species_ontologies["human"]
            sex_marker_gene_ontology_id = fixtures.validation_gene_ontologies["XIST"]
            female_ontology_id = fixtures.validation_sex_ontologies["female"]
            male_ontology_id = fixtures.validation_sex_ontologies["male"]
            MALAT1_ont_id = fixtures.validation_gene_ontologies["MALAT1"]
            human_malat1_cube = cube.df[MALAT1_ont_id, :, human_ontology_id]
            # slice cube by dimensions             gene_ontology      organ (all)          species
            human_XIST_cube = cube.df[sex_marker_gene_ontology_id, :, human_ontology_id]

            female_xist_cube = human_XIST_cube.query(f"sex_ontology_term_id == '{female_ontology_id}'")
            male_xist_cube = human_XIST_cube.query(f"sex_ontology_term_id == '{male_ontology_id}'")

            female_malat1_cube = human_malat1_cube.query(f"sex_ontology_term_id == '{female_ontology_id}'")
            male_malat1_cube = human_malat1_cube.query(f"sex_ontology_term_id == '{male_ontology_id}'")

            # should be expressed in most female cells and no male cells
            if male_xist_cube.nnz.sum() > female_xist_cube.nnz.sum():
                self.errors.append(
                    "The number of male cells expressing XIST is higher than the number of female "
                    "cells expressing XIST"
                )

            # should be expressed in females at a much higher rate. To ensure an accurate comparison divide
            # the xist expression level by the number of cells of the correct sex expressing a highly expressed
            # housekeeping gene (MALAT1 here)
            female_avg_xist_expression = female_xist_cube["sum"].sum() / female_malat1_cube["nnz"].sum()
            male_avg_xist_expression = male_xist_cube["sum"].sum() / male_malat1_cube["nnz"].sum()
            logger.info(f"female avg xist expression {female_avg_xist_expression}")
            logger.info(f"male avg xist expression {male_avg_xist_expression}")
            if male_avg_xist_expression * 50 > female_avg_xist_expression:
                self.errors.append("XIST levels dont show expected sex based difference")

    def validate_lung_cell_marker_genes(self):
        """
        Cells taken from human lungs have marker genes that are highly expressed in certain cell types
        gene: cell_type
        FCN1: monocytes
        TUBB4B: ciliated cells
        CD68: macrophages (alveioler)
        AQP5: goblet cells and secreting cells
        """
        human_ont_id = fixtures.validation_species_ontologies["human"]
        lung_ont_id = fixtures.validation_tissues_with_many_cell_types["lung"]

        # get avg expression value of gene for the celltype. That average should be greater than the avg for all
        # other cell types
        with tiledb.open(self.expression_summary_path) as cube:
            FCN1_ont_id = fixtures.validation_gene_ontologies["FCN1"]
            FCN1_human_lung_cube = cube.df[FCN1_ont_id, lung_ont_id, human_ont_id]
            self.validate_FCN1(FCN1_human_lung_cube)

            TUBB4B_ont_id = fixtures.validation_gene_ontologies["TUBB4B"]
            TUBB4B_human_lung = cube.df[TUBB4B_ont_id, lung_ont_id, human_ont_id]
            self.validate_TUBB4B(TUBB4B_human_lung)

            CD68_ont_id = fixtures.validation_gene_ontologies["CD68"]
            CD68_human_lung = cube.df[CD68_ont_id, lung_ont_id, human_ont_id]
            self.validate_CD68(CD68_human_lung)

            AQP5_ont_id = fixtures.validation_gene_ontologies["AQP5"]
            AQP5_human_lung = cube.df[AQP5_ont_id, lung_ont_id, human_ont_id]
            self.validate_AQP5(AQP5_human_lung)

    def validate_FCN1(self, FCN1_human_lung_cube):
        intermediate_monocyte_ontology_id = fixtures.validation_cell_types["intermediate monocytes"]
        non_classical_monocyte_ontology_id = fixtures.validation_cell_types["Non-classical monocyte"]
        classical_monocyte_ontology_id = fixtures.validation_cell_types["classical monocytes"]
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
        if FCN1_non_high_expression_avg > FCN1_high_expression_avg:
            self.errors.append("FCN1 expression levels are off")

    def validate_TUBB4B(self, TUBB4B_human_lung_cube):
        ciliated_cells_ontology_id = fixtures.validation_cell_types["ciliated cells"]
        lung_ciliated_ontology_id = fixtures.validation_cell_types["lung ciliated cell"]
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
        if TUBB4B_non_high_expression_avg > TUBB4B_high_expression_avg:
            self.errors.append("TUBB4B expression levels are off")

    def validate_CD68(self, CD68_human_lung_cube):
        macrophage_ontology_id = fixtures.validation_cell_types["macrophage"]
        alvelolar_macrophage_ont_id = fixtures.validation_cell_types["alveolar macrophage"]
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
        if CD68_non_high_expression_avg > CD68_high_expression_avg:
            self.errors.append("TUBB4B expression levels are off")

    def validate_AQP5(self, AQP5_human_lung_cube):
        goblet_cell_ontology_id = fixtures.validation_cell_types["goblet cells"]
        lung_goblet_ontology_id = fixtures.validation_cell_types["lung goblet cell"]
        respiratory_goblet_ontology_id = fixtures.validation_cell_types["respiratory goblet cell"]
        goblet_cell_type_ids = [goblet_cell_ontology_id, lung_goblet_ontology_id, respiratory_goblet_ontology_id]

        secretory_ontology_id = fixtures.validation_cell_types["secretory cells"]
        mucus_secretory_ontology_id = fixtures.validation_cell_types["mucus secreting cell"]
        bronchus_serous_ontology_id = fixtures.validation_cell_types["serous cell of epithelium of bronchus"]
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
        if AQP5_non_high_expression_avg > AQP5_high_expression_avg:
            self.errors.append("AQP5 expression levels are off")

    def validate_expression_levels_for_particular_gene_dataset(self):
        human_ont_id = fixtures.validation_species_ontologies["human"]
        human_lung_int = fixtures.validation_tissues_with_many_cell_types["lung"]
        MALAT1_ont_id = fixtures.validation_gene_ontologies["MALAT1"]
        CCL5_ont_id = fixtures.validation_gene_ontologies["CCL5"]
        with tiledb.open(self.expression_summary_path) as cube:
            MALAT1_human_lung_cube = cube.df[MALAT1_ont_id, human_lung_int, human_ont_id]
            CCL5_human_lung_cube = cube.df[CCL5_ont_id, human_lung_int, human_ont_id]

            MALAT1_expression = MALAT1_human_lung_cube.query(f"dataset_id == '{self.validation_dataset_id}'")
            CCL5_expression = CCL5_human_lung_cube.query(f"dataset_id == '{self.validation_dataset_id}'")

            malat1_expression_sum_by_cell_type = MALAT1_expression.groupby("cell_type_ontology_term_id").sum(
                numeric_only=True
            )["sum"]
            ccl5_expression_sum_by_cell_type = CCL5_expression.groupby("cell_type_ontology_term_id").sum(
                numeric_only=True
            )["sum"]

            expected_values = anndata.read_h5ad(
                f"{pathlib.Path(__file__).parent.resolve()}/3_0_0_lung_map_3de0ad6d-4378-4f62-b37b-ec0b75a50d94.h5ad"
            )
            malat_expected = expected_values.obs.assign(MALAT1=expected_values.layers["rankit"].toarray()[:, 0])
            ccl5_expected = expected_values.obs.assign(CCL5=expected_values.layers["rankit"].toarray()[:, 1])

            expected_malat1_by_cell_type = (
                malat_expected.groupby("cell_type_ontology_term_id").sum(numeric_only=True).MALAT1
            )
            expected_ccl5_by_cell_type = ccl5_expected.groupby("cell_type_ontology_term_id").sum(numeric_only=True).CCL5
            # drop ccl5 cell types with expression value of zero (to match pipeline processing)
            expected_ccl5_by_cell_type = expected_ccl5_by_cell_type[expected_ccl5_by_cell_type != 0]

            # ensure that both series have the same exact index ordering
            expected_malat1_by_cell_type = expected_malat1_by_cell_type[malat1_expression_sum_by_cell_type.index]
            expected_ccl5_by_cell_type = expected_ccl5_by_cell_type[ccl5_expression_sum_by_cell_type.index]

            malat1_comparison = expected_malat1_by_cell_type.compare(malat1_expression_sum_by_cell_type)
            ccl5_comparison = expected_ccl5_by_cell_type.compare(ccl5_expression_sum_by_cell_type)
            logger.info(malat1_comparison)
            logger.info(ccl5_comparison)

            """
            Because the expected values are computed using a slightly different formula they should be very close
            but not identical to the values produced by the pipeline (off by a .01 or less).
            Here we take the absolute value of the sum of the difference for each cell type. That number should be
            very small (less than 1).
            """
            malat1_diff_total = sum(abs(malat1_comparison.self - malat1_comparison.other))
            ccl5_diff_total = sum(abs(ccl5_comparison.self - ccl5_comparison.other))

            if malat1_diff_total > 1:
                self.errors.append(
                    f"MALAT1 expression values for dataset {self.validation_dataset_id} are further "
                    f"from expected values than they should be. Abs sum of difference is {malat1_diff_total}"
                )
            if ccl5_diff_total > 1:
                self.errors.append(
                    f"CCL5 expression values for dataset {self.validation_dataset_id} are further "
                    f"from expected values than they should be. Abs sum of difference is {ccl5_diff_total}"
                )

    def validate_dataset_counts(self):
        # todo check # of datasets in dataset folder and number from relational db
        with tiledb.open(self.cell_count_path) as cube:
            datasets = cube.df[:].dataset_id.drop_duplicates()
            dataset_count = len(datasets)
            if dataset_count < self.MIN_DATASET_COUNT:
                self.errors.append(
                    f"Not enough datasets in the cube, found {dataset_count} but we need {self.MIN_DATASET_COUNT}"
                )
            logger.info(f"{dataset_count} datasets included in {self.expression_summary_path}")
            logger.info(f"Included dataset ids are: {datasets}")


def query_pandas(dataframe, criteria):
    for key in criteria:
        attrs = [criteria[key]] if not isinstance(criteria[key], list) else criteria[key]
        if len(attrs) > 0:
            depluralized_key = key[:-1] if key[-1] == "s" else key
            dataframe = dataframe[dataframe[depluralized_key].isin(attrs)]
    return dataframe
