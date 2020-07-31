"""Create the 'curated' and 'remixed' datasets for the test dataset provided by Anton Ogorodnikov
from the Ye Lab. Because there is no publication or preprint yet, certain information such as
contributors and publication doi are missing. Expect an updated h5ad with an unfiltered matrix
some time in August 2020.

In terms of gene conversion, I could not merge raw counts due to the complexity of matrix
normalization. Therefore there is a duplicate gene index for the gene 'AREG' in the remixed
dataset.

The contributors described the normalization process as follows:

First ran demuxlet to remove interindividual doublets and then ran scrublet to remove
intraindividual doublets. Then we just did the following from what I remember:

####################
# Basic processing #
####################
adata = sc.read(savepath_nonorm, cache=False)
logging.info(str('Data structure details: ' + str(adata)))
# Extract list of genes
genelist = adata.var_names.tolist()
# Find mitochondrial genes
mito_genes_names = [gn for gn in genelist if gn.startswith('MT-')]
logging.info(str('Mito genes: ' + str(mito_genes_names)))
# Find indices of mitochondrial genes
mito_genes = [genelist.index(gn) for gn in mito_genes_names]
# For each cell compute fraction of counts in mito genes vs. all genes
adata.obs['percent_mito'] = np.ravel(np.sum(adata[:, mito_genes].X, axis=1)) / np.ravel(np.sum(adata.X, axis=1))
# Add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = np.ravel(adata.X.sum(axis=1))
logging.info('Filtering cells')
# Filter cells that have more than 10% of counts coming from mitochondrial genes.
adata = adata[adata.obs['percent_mito'] < 0.10]
logging.info(str('Data structure details: ' + str(adata)))
logging.info('Saving raw counts')
logging.info('Saving log(counts)+1 in .raw')
adata.raw = adata
# Filter cells with abnormally low gene counts, high gene counts.
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_cells(adata, max_genes=4000)
sc.pp.filter_genes(adata, min_cells=int(np.round(len(adata)*0.03)))
logging.info('Normalizing total counts')
sc.pp.normalize_total(adata, exclude_highly_expressed=True)
logging.info('Log transforming data')
sc.pp.log1p(adata)
logging.info(str('Data structure details: ' + str(adata)))
####################
#    Run Combat    #
####################
logging.info('Running combat on batch_cov')
sc.pp.combat(adata, key="batch_cov", covariates=["disease_cov"])
# Scale data
sc.pp.scale(adata)
logging.info('Keep highly variable genes')
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]
logging.info(str('Data structure details: ' + str(adata)))

"""

import numpy as np
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

    adata.uns['contributors'] = [{'name': 'First Author'},
                                 {'name': 'Penultimate Author', 'email': 'email@university.edu'},
                                 {'name': 'Last Author', 'email': 'email@university.edu'}]

    adata.uns['preprint_doi'] = "doi:preprint"
    adata.uns['default_embedding'] = 'X_umap'

    # Add continuous metadata
    metrics = sc.pp.calculate_qc_metrics(adata.raw, inplace=False)
    ad.obs['n_counts'] = metrics[0]['total_counts']
    ad.obs['n_genes'] = metrics[0]['n_genes_by_counts']


def remix(adata, title: str):
    """Create the full Corpora remix"""

    adata.obs[
        'assay_ontology'] = "EFO:0008995"  # ambiguous 10x assay
    adata.obs["assay"] = utils.ontology.get_ontology_label(
        "EFO:0008995")  # change from author's "10x_chromium" to "10x sequencing"
    del adata.obs["EFO"]

    disease_ontology_map = {"0000000": "PATO:0000461", "0007915": "MONDO:0007915"}
    adata.obs["disease_ontology"] = adata.obs["MONDO"].map(disease_ontology_map)
    disease_map = {disease_ontology: utils.ontology.get_ontology_label(disease_ontology) for
                   disease_ontology in adata.obs['disease_ontology'].unique()}
    adata.obs["disease"] = adata.obs["disease_ontology"].map(disease_map)
    del adata.obs["MONDO"]

    adata.obs["tissue_ontology"] = "UBERON:0001969"
    adata.obs["tissue"] = utils.ontology.get_ontology_label("UBERON:0001969")
    del adata.obs["UBERON"]
    del adata.obs["tissue_type"]

    adata.uns["organism_ontology"] = "NCBITaxon:9606"
    adata.uns["organism"] = utils.ontology.get_ontology_label("NCBITaxon:9606")

    adata.uns["title"] = title

    adata.uns["project_name"] = "Lupus single cell"
    adata.uns["project_description"] = "project description"
    adata.uns["project_raw_data_links"] = []
    adata.uns["project_other_links"] = []

    # Get cell_type from CL term, complicated by fact that one CL term is "NA" for the cell type
    # Anton called "Prolif"
    cell_type_ontology_map = {cl_term: "CL:" + cl_term for cl_term in adata.obs["CL"]}
    adata.obs["cell_type_ontology"] = adata.obs["CL"].map(cell_type_ontology_map)
    cell_type_map = {
        cl_term: utils.ontology.get_ontology_label(cl_term) if cl_term != "CL:NA" else "Prolif" for
        cl_term in adata.obs["cell_type_ontology"].unique()}
    del adata.obs["cell_type"]
    adata.obs["cell_type"] = adata.obs["cell_type_ontology"].map(cell_type_map)
    del adata.obs["CL"]

    ethnicity_ontology_map = {hancestro_term: "HANCESTRO:" + hancestro_term for hancestro_term in
                              adata.obs["HANCESTRO"]}
    adata.obs["ethnicity_ontology"] = adata.obs["HANCESTRO"].map(ethnicity_ontology_map)
    del adata.obs["HANCESTRO"]
    ethnicity_map = {ethnicity_ontology: utils.ontology.get_ontology_label(ethnicity_ontology) for
                     ethnicity_ontology in adata.obs["ethnicity_ontology"].unique()}
    adata.obs["ethnicity"] = adata.obs["ethnicity_ontology"].map(ethnicity_map)
    del adata.obs["Ethnicity"]

    # optional
    # adata.uns['tags'] = []

    # Now translate the gene symbols, which luckily can all be upgraded
    upgraded_var_index = utils.hgnc.get_upgraded_var_index(adata.var)

    # Raise exception if not all genes can be upgraded 1 to 1
    if bool(set([x for x in upgraded_var_index.to_list() if upgraded_var_index.to_list().count(x) > 1])):
        raise Exception("Need to sum gene counts after gene conversion")

    adata.var.set_index(upgraded_var_index)

    # Finally describe the layers and we're done
    adata.uns["layer_descriptions"] = {
        "raw.X": "raw",
        "X": "filter; normalize_total; log1p; combat",
    }

    return adata


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
ad = sc.read_h5ad(
    "YeLab_corpora_test.h5ad"
)
basic_curation(ad)
print_summary(ad)
ad.write("YeLab_corpora_test-curated.h5ad", compression="gzip")

rad = remix(ad, title="Ye Lab Lupus Test")
print_summary(rad)
rad.write("YeLab_corpora_test-remixed.h5ad", compression="gzip")
