"""
Create the 'curated' and 'remixed' datasets for the Strand Lab Adult Human Prostate and Urethra
Atlas (normal and diseased). Sample ethnicity field is omitted because information was missing.
There is also no cell ontology term for 'Club' and 'Hillock' epithelial cells because they are
novel cell types. Finally, the data is missing a raw count matrix so genes could not be converted.

Expect an updated h5ad from Gervaise where he will include a raw count matrix and correct a minor
bug that over-filtered some cells. The ethnicities can be added at that time.

Note from Gervaise Henry:

Data was normalized and scaled using the newer SCTransform method which does not use a scaling
factor method (uses Pearson’s residuals instead):

sc10x <- SCTransform(sc10x,vars.to.regress=c("nFeature_RNA","percent.mito"),
verbose=FALSE,return.only.var.genes=FALSE,assay="RNA")

The ethnicities are:
Sample prefix D17: Caucasian
Sample prefix D27: Caucasian
Sample prefix D35: Caucasian
Sample prefix BPH327: Caucasian
Sample prefix BPH340: Asian
Sample prefix BPH342: Caucasian
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

    # Sex of patients
    adata.obs['sex'] = "male"

    disease_ontology_map = {'BPH': 'MONDO:0010811',
                            'Normal': 'PATO:0000461'}
    adata.obs["disease_ontology"] = adata.obs["Phenotype"].map(disease_ontology_map)
    disease_map = {disease_ontology: utils.ontology.get_ontology_label(disease_ontology)
                   for disease_ontology in adata.obs["disease_ontology"].unique()}
    adata.obs["disease"] = adata.obs["disease_ontology"].map(disease_map)
    del adata.obs["Phenotype"]

    adata.obs["tissue_ontology"] = "UBERON:0002367"
    adata.obs["tissue"] = utils.ontology.get_ontology_label("UBERON:0002367")

    adata.uns["organism_ontology"] = "NCBITaxon:9606"
    adata.uns["organism"] = utils.ontology.get_ontology_label("NCBITaxon:9606")

    adata.uns['title'] = title

    adata.uns["project_name"] = "Normal and BPH Human Prostate Cell Atlas (Epithelial Cells)"
    adata.uns[
        "project_description"] = "Single-cell RNA-sequencing of adult human prostates and urethra (normal and diseased)"
    adata.uns["project_raw_data_links"] = [
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145928"]
    adata.uns["project_other_links"] = ["https://strandlab.net/sc.data/"]

    # Set the cell ontology values
    {'BE', 'Club', 'Hillock', 'LE', 'NE'}

    # Set the cell ontology values
    cell_type_map = {'BE': 'basal cell of prostate epithelium',
                     'Club': 'club cell of prostate epithelium',
                     'Hillock': 'hillock cell of prostate epithelium',
                     'LE': 'luminal cell of prostate epithelium',
                     'NE': 'endocrine-paracrine cell of prostate gland'}

    cell_type_ontology_map = {'BE': 'CL:0002341',
                              'Club': 'NA',
                              'Hillock': 'NA',
                              'LE': 'CL:0002340',
                              'NE': 'CL:0002313'}

    adata.obs["cell_type"] = adata.obs["Population"].map(cell_type_map)
    adata.obs["cell_type_ontology"] = adata.obs["Population"].map(cell_type_ontology_map)
    del adata.obs["Population"]

    # optional
    adata.uns['tags'] = ['GUDMAP', 'benign prostatic hyperplasia', 'flow cytometry',
                         'human cell atlas', 'human prostate', 'prostate cancer',
                         'prostate epithelia', 'prostate stroma', 'single-cell RNA sequencing',
                         'zonal anatomy']

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
ad = sc.read_h5ad("2020Prostate_huPr_PdPgb_epi.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("2020Prostate_huPr_PdPgb_epi-curated.h5ad", compression="gzip")

rad = remix(ad, title="Normal and BPH Human Prostate Cell Atlas (Epithelial Cells)")
print_summary(rad)
rad.write("2020Prostate_huPr_PdPgb_epi-remixed.h5ad", compression="gzip")
