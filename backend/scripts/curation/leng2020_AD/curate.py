"""Create the 'original' and 'remix' datasets for the snRNAseq of human neurons AD (Leng,
et. al. 2020) biorxiv preprint submission"""


import anndata
import numpy as np
import pandas as pd
import scanpy as sc

import utils.hgnc
import utils.ontology


def basic_curation(adata):
    """Changes to create the matrix for presentation in cellxgene."""

    # Check if there are duplicate cell or gene IDs
    if not adata.obs.index.is_unique:
        raise Exception("Cell IDs not unique.")
    if not adata.var.index.is_unique:
        raise Exception("Gene symbols not unique.")

    # These are deleted at the request of the submitter
    del adata.obsm['X_CCA']
    del adata.obsm['X_CCA.ALIGNED']

    adata.uns['contributors'] = [
        {'name': 'Kun Leng'},
        {'name': 'Emmy Li'},
        {'name': 'Rana Eser'},
        {'name': 'Antonia Piergies'},
        {'name': 'Rene Sit'},
        {'name': 'Michelle Tan'},
        {'name': 'Norma Neff'},
        {'name': 'Song Hua Li'},
        {'name': 'Roberta Diehl Rodriguez'},
        {'name': 'Claudia Kimie Suemoto'},
        {'name': 'Renata Elaine Paraizo Leite'},
        {'name': 'Carlos A. Pasqualucci'},
        {'name': 'William W. Seeley'},
        {'name': 'Salvatore Spina'},
        {'name': 'Helmut Heinsen'},
        {'name': 'Lea T. Grinberg', 'email': 'lea.grinberg@ucsf.edu'},
        {'name': 'Martin Kampmann', 'email': 'martin.kampmann@ucsf.edu'}]

    adata.uns['preprint_doi'] = "https://doi.org/10.1101/2020.04.04.025825"
    adata.uns['default_embedding'] = 'X_tSNE'


def remix(adata, title: str):
    """Create the full Corpora remix"""

    # First fill in missing metadata fields
    adata.obs['assay_ontology'] = "EFO:0009899"
    adata.obs["assay"] = utils.ontology.get_ontology_label("EFO:0009899")

    adata.obs['sex'] = "male"

    adata.obs["disease_ontology"] = "MONDO:0004975"
    adata.obs["disease"] = utils.ontology.get_ontology_label("MONDO:0004975")

    adata.obs["tissue_ontology"] = "UBERON:0002728"
    adata.obs["tissue"] = utils.ontology.get_ontology_label("UBERON:0002728")

    adata.uns["organism_ontology"] = "NCBITaxon:9606"
    adata.uns["organism"] = utils.ontology.get_ontology_label("NCBITaxon:9606")

    adata.uns["title"] = title

    adata.uns[
        "project_name"] = "Molecular characterization of selectively vulnerable neurons in " \
                          "Alzheimer’s Disease"
    adata.uns[
        "project_description"] = "Single-nuclei RNA sequencing of caudal entorhinal cortex and " \
                                 "superior frontal gyrus from individuals spanning the " \
                                 "neuropathological progression of AD"
    adata.uns["project_raw_data_links"] = [
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147528"]
    adata.uns["project_other_links"] = ["https://www.synapse.org/#!Synapse:syn21788402/wiki/601825"]

    # Set the cell ontology values
    cell_type_map = {'Exc': 'excitatory neuron',
                     'OPC': 'oligodendrocyte precursor cell',
                     'Inh': 'inhibitory neuron',
                     'Micro': 'mature microglial cell',
                     'Astro': 'mature astrocyte',
                     'Oligo': 'oligodendrocyte',
                     'Endo': 'endothelial cell'}

    adata.obs["cell_type"] = adata.obs["clusterAssignment"].str.split(":|\\.", expand=True)[1].map(
        cell_type_map)
    del adata.obs["clusterAssignment"]

    # make dictionary mapping cell_type to CL term
    cell_type_ontology_map = {cell_type: utils.ontology.lookup_candidate_term(cell_type)[0][0] for
                              cell_type in adata.obs['cell_type'].unique()}
    # result: {'excitatory neuron': 'CL:0008030', 'oligodendrocyte precursor cell': 'CL:0002453',
    # 'inhibitory neuron': 'CL:0008029', 'mature microglial cell': 'CL:0002629',
    # 'mature astrocyte': 'CL:0002627', 'oligodendrocyte': 'CL:0000128', 'endothelial cell':
    # 'CL:0000115'}

    adata.obs["cell_type_ontology"] = adata.obs["cell_type"].map(cell_type_ontology_map)

    # optional
    adata.uns['tags'] = ['AD', "Alzheimer's Disease", 'neurons']

    # Now translate the gene symbols and sum new duplicates
    # Note that we're pulling from raw here. That's where the raw counts that we can sum are
    upgraded_var_index = utils.hgnc.get_upgraded_var_index(adata.var)
    merged_raw_counts = pd.DataFrame.sparse.from_spmatrix(
        adata.X, index=adata.obs.index, columns=upgraded_var_index,
    ).sum(axis=1, level=0, skipna=False)

    # Create the new anndata object with the summed values
    remix_adata = anndata.AnnData(
        X=merged_raw_counts,
        obs=adata.obs,
        var=merged_raw_counts.columns.to_frame(name="hgnc_gene_symbol"),
        uns=adata.uns,
        obsm=adata.obsm,
    )

    # Finally describe the layers and we're done
    remix_adata.uns["layer_descriptions"] = {
        "X": "raw"
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


# Process EC_all
ad = sc.read_h5ad("EC_allCells/kampmann_lab_human_AD_snRNAseq_EC.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("EC_allCells/kampmann_lab_human_AD_snRNAseq_EC-curated.h5ad", compression="gzip")
rad = remix(ad, title="Molecular characterization of selectively vulnerable neurons in "
                      "Alzheimer’s Disease: caudal entorhinal cortex")
print_summary(rad)
rad.write("EC_allCells/kampmann_lab_human_AD_snRNAseq_EC-remixed.h5ad", compression="gzip")

# Process SFG_all
ad = sc.read_h5ad("SFG_allCells/kampmann_lab_human_AD_snRNAseq_SFG.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("SFG_allCells/kampmann_lab_human_AD_snRNAseq_SFG-curated.h5ad", compression="gzip")
rad = remix(ad, title="Molecular characterization of selectively vulnerable neurons in "
                      "Alzheimer’s Disease: superior frontal gyrus")
print_summary(rad)
rad.write("SFG_allCells/kampmann_lab_human_AD_snRNAseq_SFG-remixed.h5ad", compression="gzip")

# Process EC_astrocytes
ad = sc.read_h5ad("EC_subclusters/EC_astrocytes/kampmann_lab_human_AD_snRNAseq_EC_astrocytes.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("EC_subclusters/EC_astrocytes/kampmann_lab_human_AD_snRNAseq_EC_astrocytes-curated.h5ad", compression="gzip")
rad = remix(ad, title="Molecular characterization of selectively vulnerable neurons in "
                      "Alzheimer’s Disease: EC astrocytes")
print_summary(rad)
rad.write("EC_subclusters/EC_astrocytes/kampmann_lab_human_AD_snRNAseq_EC_astrocytes-remixed.h5ad", compression="gzip")

# Process EC_excitatoryNeurons
ad = sc.read_h5ad("EC_subclusters/EC_excitatoryNeurons/kampmann_lab_human_AD_snRNAseq_EC_excitatoryNeurons.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("EC_subclusters/EC_excitatoryNeurons/kampmann_lab_human_AD_snRNAseq_EC_excitatoryNeurons-curated.h5ad",
         compression="gzip")
rad = remix(ad, title="Molecular characterization of selectively vulnerable neurons in "
                      "Alzheimer’s Disease: EC excitatoryNeurons")
print_summary(rad)
rad.write("EC_subclusters/EC_excitatoryNeurons/kampmann_lab_human_AD_snRNAseq_EC_excitatoryNeurons-remixed.h5ad",
          compression="gzip")

# Process EC_inhibitoryNeurons
ad = sc.read_h5ad("EC_subclusters/EC_inhibitoryNeurons/kampmann_lab_human_AD_snRNAseq_EC_inhibitoryNeurons.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("EC_subclusters/EC_inhibitoryNeurons/kampmann_lab_human_AD_snRNAseq_EC_inhibitoryNeurons-curated.h5ad",
         compression="gzip")
rad = remix(ad, title="Molecular characterization of selectively vulnerable neurons in "
                      "Alzheimer’s Disease: EC inhibitoryNeurons")
print_summary(rad)
rad.write("EC_subclusters/EC_inhibitoryNeurons/kampmann_lab_human_AD_snRNAseq_EC_inhibitoryNeurons-remixed.h5ad",
          compression="gzip")

# Process EC_microglia
ad = sc.read_h5ad("EC_subclusters/EC_microglia/kampmann_lab_human_AD_snRNAseq_EC_microglia.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("EC_subclusters/EC_microglia/kampmann_lab_human_AD_snRNAseq_EC_microglia-curated.h5ad", compression="gzip")
rad = remix(ad, title="Molecular characterization of selectively vulnerable neurons in "
                      "Alzheimer’s Disease: EC microglia")
print_summary(rad)
rad.write("EC_subclusters/EC_microglia/kampmann_lab_human_AD_snRNAseq_EC_microglia-remixed.h5ad", compression="gzip")

# Process SFG_astrocytes
ad = sc.read_h5ad("SFG_subclusters/SFG_astrocytes/kampmann_lab_human_AD_snRNAseq_SFG_astrocytes.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("SFG_subclusters/SFG_astrocytes/kampmann_lab_human_AD_snRNAseq_SFG_astrocytes-curated.h5ad",
         compression="gzip")
rad = remix(ad, title="MolSFGular characterization of selSFGtively vulnerable neurons in "
                      "Alzheimer’s Disease: SFG astrocytes")
print_summary(rad)
rad.write("SFG_subclusters/SFG_astrocytes/kampmann_lab_human_AD_snRNAseq_SFG_astrocytes-remixed.h5ad",
          compression="gzip")

# Process SFG_excitatoryNeurons
ad = sc.read_h5ad("SFG_subclusters/SFG_excitatoryNeurons/kampmann_lab_human_AD_snRNAseq_SFG_excitatoryNeurons.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("SFG_subclusters/SFG_excitatoryNeurons/kampmann_lab_human_AD_snRNAseq_SFG_excitatoryNeurons-curated.h5ad",
         compression="gzip")
rad = remix(ad, title="MolSFGular characterization of selSFGtively vulnerable neurons in "
                      "Alzheimer’s Disease: SFG excitatoryNeurons")
print_summary(rad)
rad.write("SFG_subclusters/SFG_excitatoryNeurons/kampmann_lab_human_AD_snRNAseq_SFG_excitatoryNeurons-remixed.h5ad",
          compression="gzip")

# Process SFG_inhibitoryNeurons
ad = sc.read_h5ad("SFG_subclusters/SFG_inhibitoryNeurons/kampmann_lab_human_AD_snRNAseq_SFG_inhibitoryNeurons.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("SFG_subclusters/SFG_inhibitoryNeurons/kampmann_lab_human_AD_snRNAseq_SFG_inhibitoryNeurons-curated.h5ad",
         compression="gzip")
rad = remix(ad, title="MolSFGular characterization of selSFGtively vulnerable neurons in "
                      "Alzheimer’s Disease: SFG inhibitoryNeurons")
print_summary(rad)
rad.write("SFG_subclusters/SFG_inhibitoryNeurons/kampmann_lab_human_AD_snRNAseq_SFG_inhibitoryNeurons-remixed.h5ad",
          compression="gzip")

# Process SFG_microglia
ad = sc.read_h5ad("SFG_subclusters/SFG_microglia/kampmann_lab_human_AD_snRNAseq_SFG_microglia.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("SFG_subclusters/SFG_microglia/kampmann_lab_human_AD_snRNAseq_SFG_microglia-curated.h5ad", compression="gzip")
rad = remix(ad, title="MolSFGular characterization of selSFGtively vulnerable neurons in "
                      "Alzheimer’s Disease: SFG microglia")
print_summary(rad)
rad.write("SFG_subclusters/SFG_microglia/kampmann_lab_human_AD_snRNAseq_SFG_microglia-remixed.h5ad",
          compression="gzip")