"""Create the 'original' and 'remix' datasets for the scv2 submission"""


import anndata
import numpy as np
import pandas as pd
import scanpy as sc

import utils.hgnc
import utils.ontology


def basic_curation(adata):
    """Changes to create the matrix for presenation in cellxgene."""

    # These are deleted at the request of the submitter
    del adata.obs["scv2+"]
    del adata.obs["scv2_5+"]
    del adata.obs["res_inf"]
    del adata.obs["ees_inf"]
    del adata.obs["subctype"]
    del adata.obs["init_ctype"]
    del adata.obs["raw"]
    del adata.obs["batch"]
    del adata.obs["ctype_res0.5"]
    del adata.obs["louvain_res0.5"]

    # This is also some logic the submitter requested
    adata.obs["Infected"] = adata.obs["scv2_10+"].astype("category").replace({1: "Infected", 0: "Bystander"})
    del adata.obs["scv2_10+"]
    adata.obs["Infected"][adata.obs["Condition"] == "Mock"] = "Mock"

    adata.obs["Cell type"] = adata.obs["ctype"]
    del adata.obs["ctype"]

    adata.uns["contributors"] = [
        {"name": "Neal G. Ravindra", "institution": "Yale University", "email": "neal.ravindra@yale.edu"}
    ]
    adata.uns["preprint_doi"] = "https://doi.org/10.1101/2020.05.06.081695"


def remix(adata):
    """Create the full Corpora remix"""

    # First fill in missing metadata fields
    adata.obs["assay"] = "10x 3' v3 sequencing"
    adata.obs["assay_ontology"] = "EFO:0009922"

    adata.obs["disease"] = adata.obs["Condition"].replace(
        {"1dpi": "SARS-CoV-2", "2dpi": "SARS-CoV-2", "3dpi": "SARS-CoV-2", "Mock": "normal"}, inplace=False,
    )
    adata.obs["disease_ontology"] = adata.obs["Condition"].replace(
        {"1dpi": "MONDO:0100096", "2dpi": "MONDO:0100096", "3dpi": "MONDO:0100096", "Mock": "PATO:0000461"},
        inplace=False,
    )
    adata.obs["tissue"] = "epithelium of bronchus (cell culture)"
    adata.obs["tissue_ontology"] = "UBERON:0002031 (cell culture)"

    adata.uns["organism"] = "Homo sapiens"
    adata.uns["organism_ontology"] = "NCBITaxon:9606"
    adata.uns["title"] = (
        "Single-cell longitudinal analysis of SARS-CoV-2 infection in human " "bronchial epithelial cells"
    )

    adata.uns["project_name"] = (
        "Single-cell longitudinal analysis of SARS-CoV-2 infection in " "human bronchial epithelial cells"
    )
    adata.uns["project_description"] = (
        "Single-cell RNA sequencing of experimentally infected human bronchial epithelial "
        "cells (HBECs) in air-liquid interface cultures over a time-course"
    )
    adata.uns["project_other_links"] = ["https://github.com/vandijklab/HBEC_SARS-CoV-2_scRNA-seq"]
    # Set the cell ontology values
    # BC/Club is "a cell population intermediate between basal cells and club cells". There is no
    # CL term for them.

    cell_type_ontology_map = {
        "Basal cells": "CL:1000349",
        "Ciliated cells": "CL:0002332",
        "Club cells": "CL:0000158",
        "BC/Club": "",
        "Ionocytes": "CL:0017000",
        "Neuroendocrine cells": "CL:1000223",
        "Tuft cells": "CL:0002208",
        "Goblet cells": "CL:1000312",
    }
    cell_type_map = {k: utils.ontology.get_ontology_label(v) for k, v in cell_type_ontology_map.items()}
    cell_type_map["BC/Club"] = "BC/Club"

    adata.obs["cell_type_ontology"] = adata.obs["Cell type"].replace(cell_type_ontology_map)
    adata.obs["cell_type"] = adata.obs["Cell type"].replace(cell_type_map)
    del adata.obs["Cell type"]

    # Now translate the gene symbols and sum new duplicates
    # Note that we're pulling from raw here. That's where the raw counts that we can sum are
    upgraded_var_index = utils.hgnc.get_upgraded_var_index(adata.var)
    merged_raw_counts = pd.DataFrame.sparse.from_spmatrix(
        adata.raw.X, index=adata.obs.index, columns=upgraded_var_index,
    ).sum(axis=1, level=0, skipna=False)

    # Create the new anndata object with the summed values
    remix_adata = anndata.AnnData(
        X=merged_raw_counts,
        obs=adata.obs,
        var=merged_raw_counts.columns.to_frame(name="hgnc_gene_symbol"),
        uns=adata.uns,
        obsm=adata.obsm,
    )
    remix_adata.raw = remix_adata.copy()

    # Perform the same tranformations on the new values as they did in the paper
    sc.pp.normalize_total(remix_adata)
    sc.pp.sqrt(remix_adata)

    # Finally describe the layers and we're done
    remix_adata.uns["layer_descriptions"] = {
        "raw.X": "raw",
        "X": "normalize_total; sqrt",
    }

    return remix_adata


def print_summary(adata):
    """Print out a little summary of the metadata."""
    print(adata.obs.dtypes)
    for column in adata.obs.nunique().items():
        field, n_unique = column
        if n_unique > 1000 and not np.issubdtype(adata.obs[field].dtype, np.number):
            print("TOO MANY:", field)
        elif n_unique == 1:
            print("ONLY ONE:", field)


ad = sc.read_h5ad(
    "Single_cell_longitudinal_analysis_of_SARS_CoV_2_infection_in_human_bronchial_epithelial_cells-29-original.h5ad"
)
basic_curation(ad)
print_summary(ad)
ad.write(
    "Single_cell_longitudinal_analysis_of_SARS_CoV_2_infection_in_human_bronchial_epithelial_cells-29-curated.h5ad",
    compression="gzip",
)
rad = remix(ad)
print_summary(rad)
rad.write(
    "Single_cell_longitudinal_analysis_of_SARS_CoV_2_infection_in_human_bronchial_epithelial_cells-29-remixed.h5ad",
    compression="gzip",
)
