import logging

import anndata
from pathlib import Path

def validate_corpus_load(anndata_object: anndata.AnnData, group_name: str, dataset_id: str):
    """
    Validate that the load looks sane
    """
    pass


def validate_cube(snapshot):
    pass

# grab cube
# check size
# check size of human v mouse
# list size of tissues?

# check MALAT1 and ACTB

# check XIST appears in women but not men

# check human lung cells of particular types have marker genes
# FCN1: monocytes
# TUBB4B: ciliated cells, and this gene is expressed to lower values in most other cell types
# CD68: macrophages
# AQP5: goblet cells and secreting cells

# check expression levels are correct for lung map dataset uuid 3de0ad6-4378-4f62-b37b-ec0b75a50d94
# genes ["MALAT1", "CCL5"]

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


validation_ontologies = {

}

def validate_cube_size(path_to_cube):
    size_byte=sum(file.stat().st_size for file in Path('expression_summary').rglob('*'))
    size_gb = size_byte/1000/1024/1024
    # check expression cube is larger than 3gb
    assert(size_gb>3)
    logger.info(f"Expression summary cube is {size_gb}GB")


def validate_cube_species():
    pass


def validate_tissues_in_cube():
    pass


def validate_housekeeping_gene_expression_levels():
    pass


def validate_sex_specific_marker_gene():
    pass


def validate_lung_cell_marker_genes():
    pass


def validate_expression_levels_for_particular_gene_dataset():
    pass
