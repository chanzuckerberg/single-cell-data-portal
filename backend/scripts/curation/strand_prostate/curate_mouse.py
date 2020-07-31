"""
Create the 'curated' and 'remixed' datasets for the Strand Lab Mouse Prostate and Urethral Cell
Atlas. Sample ethnicity field is omitted because mouse. There was also no cell ontology term for
'urethral luminal' cells. Finally, the data is missing a raw count matrix so genes could not be
converted even if we had a standard set of mouse genes.

Note from Gervaise Henry:

Data was normalized and scaled using the newer SCTransform method which does not use a scaling
factor method (uses Pearson’s residuals instead):

sc10x <- SCTransform(sc10x,vars.to.regress=c("nFeature_RNA","percent.mito"),
verbose=FALSE,return.only.var.genes=FALSE,assay="RNA")

Expect an updated h5ad from Gervaise where he will include a raw count matrix and correct a minor
bug that over-filtered some cells.
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

    adata.uns['contributors'] = [
        {'name': 'Diya B. Joseph'},
        {'name': 'Gervaise H. Henry'},
        {'name': 'Alicia Malewska'},
        {'name': 'Nida S. Iqbal'},
        {'name': 'Hannah M. Ruetten'},
        {'name': 'Anne E. Turco'},
        {'name': 'Lisa L. Abler'},
        {'name': 'Simran K. Sandhu'},
        {'name': 'Mark T. Cadena'},
        {'name': 'Venkat S. Malladi'},
        {'name': 'Jeffrey C. Reese'},
        {'name': 'Ryan J. Mauck'},
        {'name': 'Jeffrey C. Gahan'},
        {'name': 'Ryan C. Hutchinson'},
        {'name': 'Claus G. Roehrborn'},
        {'name': 'Linda A. Baker'},
        {'name': 'Chad M. Vezina'},
        {'name': 'Douglas W Strand', 'email': 'douglas.strand@utsouthwestern.edu'}]

    adata.uns['preprint_doi'] = "https://doi.org/10.1101/2020.02.19.937615"
    adata.uns['publication_doi'] = "https://doi.org/10.1002/pros.24020"

    adata.uns['default_embedding'] = 'X_umap'


def remix(adata, title: str):
    """Create the full Corpora remix"""

    # First fill in missing metadata fields
    adata.obs['assay_ontology'] = "EFO:0009899"
    adata.obs["assay"] = utils.ontology.get_ontology_label("EFO:0009899")

    adata.obs['sex'] = "male"

    adata.obs["disease_ontology"] = "PATO:0000461"
    adata.obs["disease"] = utils.ontology.get_ontology_label("PATO:0000461")

    tissue_ontology_map = {'Prostate': 'UBERON:0002367',
                           'Urethra': 'UBERON:0001333'}
    adata.obs['tissue_ontology'] = adata.obs['Region'].map(tissue_ontology_map)
    del adata.obs['Region']

    tissue_map = {uberon: utils.ontology.get_ontology_label(uberon) for
                  uberon in ad.obs['tissue_ontology'].unique()}
    adata.obs['tissue'] = adata.obs['tissue_ontology'].map(tissue_map)

    adata.uns["organism_ontology"] = "NCBITaxon:10090"
    adata.uns["organism"] = utils.ontology.get_ontology_label("NCBITaxon:10090")

    adata.uns['title'] = title

    adata.uns["project_name"] = "Mouse Prostate and Urethral Cell Atlas"
    adata.uns[
        "project_description"] = "Single-cell RNA-sequencing of adult mouse lower urinary tracts."
    adata.uns["project_raw_data_links"] = [
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145929"]
    adata.uns["project_other_links"] = [
        "https://www.gudmap.org/chaise/record/#2/Common:Collection/RID=16-WM8C",
        "https://strandlab.net/sc.data/"]

    # Set the cell ontology values
    cell_type_map = {'BE': 'basal cell of prostate epithelium',
                     'LE': 'luminal cell of prostate epithelium',
                     'Ur': 'urethral luminal'}

    cell_type_ontology_map = {'BE': 'CL:0002341',
                              'LE': 'CL:0002340',
                              'Ur': 'NA'}

    adata.obs["cell_type"] = adata.obs["Population"].map(cell_type_map)
    adata.obs["cell_type_ontology"] = adata.obs["Population"].map(cell_type_ontology_map)
    del adata.obs["Population"]

    # optional
    adata.uns['tags'] = ['GUDMAP', 'flow cytometry', 'mouse cell atlas', 'mouse prostate',
                         'prostate cancer', 'prostate epithelia', 'prostate stroma',
                         'single-cell RNA sequencing', 'zonal anatomy']

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


# Process mouse epithelial h5ad
ad = sc.read_h5ad("2020Prostate_muPrUr_epi.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("2020Prostate_muPrUr_epi-curated.h5ad", compression="gzip")

rad = remix(ad, title="Mouse Prostate and Urethral Cell Atlas (Epithelial Cells)")
print_summary(rad)
rad.write("2020Prostate_muPrUr_epi-remixed.h5ad", compression="gzip")
