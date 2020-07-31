"""Create the 'remix' dataset for the Krasnow Lab Human Lung Cell Atlas, 10x. The cell type and
cell ontology was manually curated and can be seen in the cell_type_map and
cell_type_ontology_map dicts."""


import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import json

import utils.hgnc
import utils.ontology


def basic_curation(adata):
    """Changes to create the matrix for presentation in cellxgene."""

    # Check if there are duplicate cell or gene IDs
    if not adata.obs.index.is_unique:
        raise Exception("Cell IDs not unique.")
    if not adata.var.index.is_unique:
        raise Exception("Gene symbols not unique.")

    # Contributors
    # must be list of dictionaries but original data 'authors' is list of string represetations
    # of dictionaries
    # first make strings json acceptable while maintaining compatibility with names that contain
    # apostrophe
    json_acceptable_str = [author_dict.replace("{\'", "{\"").replace("\': \'", "\": \"").replace(
        "\', \'", "\", \"").replace("\'}", "\"}") for author_dict in list(adata.uns['authors'])]
    # decode str into dict
    adata.uns['contributors'] = [json.loads(author_dict) for author_dict in json_acceptable_str]
    del adata.uns['authors']

    adata.uns['preprint_doi'] = adata.uns['preprint']['doi']
    del adata.uns['preprint']['doi']

    adata.uns['default_embedding'] = 'X_tSNE'


def remix(adata, title: str):
    """Create the full Corpora remix"""

    # First fill in missing metadata fields
    adata.obs['assay_ontology'] = "EFO:0009899"
    adata.obs["assay"] = utils.ontology.get_ontology_label("EFO:0009899")

    # Sex of patients
    sex_dict = {1: "male", 2: "male", 3: "female"}  # map patient number to sex
    adata.obs['sex'] = adata.obs['patient'].map(sex_dict)

    adata.obs["disease_ontology"] = "PATO:0000461"
    adata.obs["disease"] = utils.ontology.get_ontology_label("PATO:0000461")

    tissue_ontology_dict = {'lung': "UBERON:0002048", 'blood': "UBERON:0000178"}
    adata.obs['tissue_ontology'] = adata.obs['tissue'].map(tissue_ontology_dict)

    adata.uns["organism_ontology"] = "NCBITaxon:9606"
    adata.uns["organism"] = utils.ontology.get_ontology_label("NCBITaxon:9606")

    # title
    adata.uns['title'] = title
    del adata.uns['preprint']['title']

    adata.uns[
        "project_name"] = "A molecular cell atlas of the human lung from single cell RNA sequencing"
    adata.uns[
        "project_description"] = "Droplet- and plate-based single cell RNA sequencing applied to " \
                                 "∼75,000 human lung and blood cells, combined with a " \
                                 "multi-pronged cell annotation approach, which have allowed us " \
                                 "to define the gene expression profiles and anatomical locations " \
                                 "of 58 cell populations in the human lung, including 41 of 45 " \
                                 "previously known cell types or subtypes and 14 new ones."
    adata.uns["project_raw_data_links"] = ["https://www.ebi.ac.uk/ega/studies/EGAS00001004344"]
    adata.uns["project_other_links"] = ["https://hlca.ds.czbiohub.org/",
                                        "https://www.synapse.org/#!Synapse:syn21041850/wiki/600865",
                                        "https://github.com/krasnowlab/HLCA"]

    # Set the cell ontology values
    cell_type_map = {'Capillary Aerocyte': 'Capillary Aerocyte',
                     'Capillary': 'Capillary Endothelial Cell',
                     'Capillary Intermediate 1': 'Capillary Intermediate 1',
                     'Capillary Intermediate 2': 'Capillary Intermediate 2',
                     'IGSF21+ Dendritic': 'IGSF21+ Dendritic',
                     'Myeloid Dendritic Type 1': 'Myeloid Dendritic Cell, Human',
                     'Plasmacytoid Dendritic': 'Plasmacytoid Dendritic Cell, Human',
                     'Myeloid Dendritic Type 2': 'CD1c-positive Myeloid Dendritic Cell',
                     'B': 'B Cell',
                     'EREG+ Dendritic': 'EREG+ Dendritic',
                     'Macrophage': 'Alveolar Macrophage',
                     'CD8+ Naive T': 'Naive Thymus-derived CD8-positive, alpha-beta T Cell',
                     'CD4+ Naive T': 'Naive Thymus-derived CD4-positive, alpha-beta T Cell',
                     'CD4+ Memory/Effector T': 'Effector Memory CD4-positive, alpha-beta T Cell',
                     'Vein': 'Vein Endothelial Cell',
                     'Artery': 'Endothelial Cell of Artery',
                     'Pericyte': 'Pericyte Cell',
                     'Vascular Smooth Muscle': 'Vascular Associated Smooth Muscle Cell',
                     'Club': 'Club Cell',
                     'Mucous': 'Mucus Secreting Cell',
                     'Alveolar Epithelial Type 2': 'Type II Pneumocyte',
                     'Basal': 'Respiratory Basal Cell',
                     'Lymphatic': 'Endothelial Cell of Lymphatic Vessel',
                     'Proliferating Macrophage': 'Proliferating Macrophage',
                     'CD8+ Memory/Effector T': 'Effector Memory CD8-positive, alpha-beta T Cell',
                     'Proliferating NK/T': 'Proliferating NK/T',
                     'Natural Killer T': 'Mature NK T Cell',
                     'Natural Killer': 'Natural Killer Cell',
                     'OLR1+ Classical Monocyte': 'OLR1+ Classical Monocyte',
                     'Basophil/Mast 1': 'Basophil/Mast 1',
                     'Classical Monocyte': 'Classical Monocyte',
                     'Intermediate Monocyte': 'Intermediate Monocyte',
                     'Nonclassical Monocyte': 'Non-classical Monocyte',
                     'Airway Smooth Muscle': 'Bronchial Smooth Muscle Cell',
                     'Ciliated': 'Lung Ciliated Cell',
                     'Alveolar Fibroblast': 'Alveolar Fibroblast',
                     'Myofibroblast': 'Myofibroblast Cell',
                     'Adventitial Fibroblast': 'Adventitial Fibroblast',
                     'Alveolar Epithelial Type 1': 'Type I Pneumocyte',
                     'Platelet/Megakaryocyte': 'Megakaryocyte',
                     'TREM2+ Dendritic': 'TREM2+ Dendritic',
                     'Differentiating Basal': 'Differentiating Basal',
                     'Proliferating Basal': 'Proliferating Basal',
                     'Plasma': 'Plasma Cell',
                     'Bronchial Vessel 2': 'Bronchial Vessel 2',
                     'Bronchial Vessel 1': 'Bronchial Vessel 1',
                     'Lipofibroblast': 'Lipofibroblast',
                     'Mesothelial': 'Mesothelial Cell of Pleura',
                     'Basophil/Mast 2': 'Basophil/Mast 2',
                     'Signaling Alveolar Epithelial Type 2': 'Signaling Alveolar Epithelial Type 2',
                     'Proximal Basal': 'Respiratory Basal Cell',
                     'Neuroendocrine': 'Lung Neuroendocrine Cell',
                     'Fibromyocyte': 'Fibromyocyte',
                     'Ionocyte': 'Pulmonary Ionocyte',
                     'Serous': 'Tracheobronchial Serous Cell',
                     'Proximal Ciliated': 'Proximal Ciliated',
                     'Goblet': 'Lung Goblet Cell'}

    adata.obs["cell_type"] = adata.obs["free_annotation"].str.split("_", expand=True)[0].map(
        cell_type_map)

    # dictionary mapping clusterAssignment to CL term
    cell_type_ontology_map = {'Capillary Aerocyte': '',
                              'Capillary': 'CL:0002144',
                              'Capillary Intermediate 1': '',
                              'Capillary Intermediate 2': '',
                              'IGSF21+ Dendritic': '',
                              'Myeloid Dendritic Type 1': 'CL:0001057',
                              'Plasmacytoid Dendritic': 'CL:0001058',
                              'Myeloid Dendritic Type 2': 'CL:0002399',
                              'B': 'CL:0000236',
                              'EREG+ Dendritic': '',
                              'Macrophage': 'CL:0000583',
                              'CD8+ Naive T': 'CL:0000900',
                              'CD4+ Naive T': 'CL:0000895',
                              'CD4+ Memory/Effector T': 'CL:0000905',
                              'Vein': 'CL:0002543',
                              'Artery': 'CL:1000413',
                              'Pericyte': 'CL:0000669',
                              'Vascular Smooth Muscle': 'CL:0000359',
                              'Club': 'CL:0000158',
                              'Mucous': 'CL:0000319',
                              'Alveolar Epithelial Type 2': 'CL:0002063',
                              'Basal': 'CL:0002633',
                              'Lymphatic': 'CL:0002138',
                              'Proliferating Macrophage': '',
                              'CD8+ Memory/Effector T': 'CL:0000913',
                              'Proliferating NK/T': '',
                              'Natural Killer T': 'CL:0000814',
                              'Natural Killer': 'CL:0000623',
                              'OLR1+ Classical Monocyte': '',
                              'Basophil/Mast 1': '',
                              'Classical Monocyte': 'CL:0000860',
                              'Intermediate Monocyte': 'CL:0002393',
                              'Nonclassical Monocyte': 'CL:0000875',
                              'Airway Smooth Muscle': 'CL:0002598',
                              'Ciliated': 'CL:1000271',
                              'Alveolar Fibroblast': '',
                              'Myofibroblast': 'CL:0000186',
                              'Adventitial Fibroblast': '',
                              'Alveolar Epithelial Type 1': 'CL:0002062',
                              'Platelet/Megakaryocyte': 'CL:0000556',
                              'TREM2+ Dendritic': '',
                              'Differentiating Basal': '',
                              'Proliferating Basal': '',
                              'Plasma': 'CL:0000786',
                              'Bronchial Vessel 2': '',
                              'Bronchial Vessel 1': '',
                              'Lipofibroblast': '',
                              'Mesothelial': 'CL:1000491',
                              'Basophil/Mast 2': '',
                              'Signaling Alveolar Epithelial Type 2': '',
                              'Proximal Basal': 'CL:0002633',
                              'Neuroendocrine': 'CL:1000223',
                              'Fibromyocyte': '',
                              'Ionocyte': 'CL:0017000',
                              'Serous': 'CL:0019001',
                              'Proximal Ciliated': '',
                              'Goblet': 'CL:1000143'}

    adata.obs["cell_type_ontology"] = adata.obs["free_annotation"].str.split("_", expand=True)[
        0].map(cell_type_ontology_map)
    del adata.obs["free_annotation"]

    # optional
    adata.uns['tags'] = ['lung', 'human cell atlas', 'HCA', 'HLCA']

    # Now translate the gene symbols and sum new duplicates
    # Note that we're pulling from raw here. That's where the raw counts that we can sum are
    upgraded_var_index = utils.hgnc.get_upgraded_var_index(adata.var)
    merged_raw_counts = pd.DataFrame(adata.raw.X, index=adata.obs.index,
                                     columns=upgraded_var_index).sum(axis=1, level=0, skipna=False)

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
    sc.pp.normalize_total(remix_adata, target_sum=1e4)
    sc.pp.log1p(remix_adata)

    # Finally describe the layers and we're done
    remix_adata.uns["layer_descriptions"] = {
        "raw.X": "raw",
        "X": "normalize_total; log1p",
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

    # Print missing cell fields required by Corpora schema
    remix_cellfields = np.array(
        ['tissue', 'assay', 'disease', 'cell_type', 'sex', 'ethnicity', 'tissue_ontology',
         'assay_ontology', 'disease_ontology', 'cell_type_ontology', 'ethnicity_ontology'])
    missing_remix_cellfields = np.array(set(remix_cellfields) - set(adata.obs.columns.values))
    print("MISSING CORPORA FIELDS:", missing_remix_cellfields)


# Process h5ad
ad = sc.read_h5ad("krasnow_lab_human_lung_cell_atlas_10x-1-curated.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("krasnow_lab_human_lung_cell_atlas_10x-1-curatedv2.h5ad", compression="gzip")
rad = remix(ad, title="Krasnow Lab Human Lung Cell Atlas, 10X")
print_summary(rad)
rad.write("krasnow_lab_human_lung_cell_atlas_10x-1-remixed.h5ad", compression="gzip")
