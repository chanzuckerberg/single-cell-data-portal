"""
Create the 'curated' and 'remixed' datasets for the Strand Lab Normal Human Prostate Cell Atlas.
Sample ethnicity field is omitted because information was missing. There is also no cell ontology
term for 'Club' and 'Hillock' epithelial cells because they are novel cell types. Finally,
the data is missing a raw count matrix so genes could not be converted.

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

    adata.uns["contributors"] = [
        {"name": "Gervaise H Henry"},
        {"name": "Alicia Malewska"},
        {"name": "Diya B Joseph"},
        {"name": "Venkat S Malladi"},
        {"name": "Jeon Lee"},
        {"name": "Jose Torrealba"},
        {"name": "Ryan J Mauck"},
        {"name": "Jeffrey C Gahan"},
        {"name": "Ganesh V Raj"},
        {"name": "Claus G Roehrborn"},
        {"name": "Gary C Hon"},
        {"name": "Malcolm P MacConmara"},
        {"name": "Jeffrey C Reese"},
        {"name": "Ryan C Hutchinson"},
        {"name": "Chad M Vezina"},
        {"name": "Douglas W Strand", "email": "douglas.strand@utsouthwestern.edu"},
    ]

    adata.uns["preprint_doi"] = "https://doi.org/10.1101/439935"
    adata.uns["publication_doi"] = "https://doi.org/10.1016/j.celrep.2018.11.086"

    adata.uns["default_embedding"] = "X_tsne"


def remix(adata, title: str):
    """Create the full Corpora remix"""

    # First fill in missing metadata fields
    adata.obs["assay_ontology"] = "EFO:0009899"
    adata.obs["assay"] = utils.ontology.get_ontology_label("EFO:0009899")

    # Sex of patients
    adata.obs["sex"] = "male"

    adata.obs["disease_ontology"] = "PATO:0000461"
    adata.obs["disease"] = utils.ontology.get_ontology_label("PATO:0000461")

    adata.obs["tissue_ontology"] = "UBERON:0002367"
    adata.obs["tissue"] = utils.ontology.get_ontology_label("UBERON:0002367")

    adata.uns["organism_ontology"] = "NCBITaxon:9606"
    adata.uns["organism"] = utils.ontology.get_ontology_label("NCBITaxon:9606")

    adata.uns["title"] = title

    adata.uns["project_name"] = "Normal Human Prostate Cell Atlas"
    adata.uns["project_description"] = (
        "Creation of a cellular anatomy of the young human prostate by scRNA " "sequencing."
    )
    adata.uns["project_raw_data_links"] = ["https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117403"]
    adata.uns["project_other_links"] = ["https://git.biohpc.swmed.edu/StrandLab/sc-TissueMapper_Pr"]

    # Set the cell ontology values
    # BE = basal epithelial
    # LE = luminal epithelial
    # NE = neuroendocrine epithelial
    # Club = club epithelial
    # Endo = endothelial
    # Fib = fibroblast
    # Hillock = hillock epithelial
    # Leu = leukocyte
    # SM = smooth muscle

    cell_type_map = {
        "BE": "basal cell of prostate epithelium",
        "Club": "club cell of prostate epithelium",
        "Endo": "prostate gland microvascular endothelial cell",
        "Fib": "fibroblast of connective tissue of prostate",
        "Hillock": "hillock cell of prostate epithelium",
        "LE": "luminal cell of prostate epithelium",
        "Leu": "leukocyte",
        "NE": "endocrine-paracrine cell of prostate gland",
        "SM": "smooth muscle cell of prostate",
    }

    cell_type_ontology_map = {
        "BE": "CL:0002341",
        "Club": "NA",
        "Endo": "CL:2000059",
        "Fib": "CL:1000299",
        "Hillock": "NA",
        "LE": "CL:0002340",
        "Leu": "CL:0000738",
        "NE": "CL:0002313",
        "SM": "CL:1000487",
    }

    adata.obs["cell_type"] = adata.obs["Population"].map(cell_type_map)
    adata.obs["cell_type_ontology"] = adata.obs["Population"].map(cell_type_ontology_map)
    del adata.obs["Population"]

    # optional
    adata.uns["tags"] = [
        "GUDMAP",
        "benign prostatic hyperplasia",
        "flow cytometry",
        "human cell atlas",
        "human prostate",
        "prostate cancer",
        "prostate epithelia",
        "prostate stroma",
        "single-cell RNA sequencing",
        "zonal anatomy",
    ]

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
        [
            "tissue",
            "assay",
            "disease",
            "cell_type",
            "sex",
            "ethnicity",
            "tissue_ontology",
            "assay_ontology",
            "disease_ontology",
            "cell_type_ontology",
            "ethnicity_ontology",
        ]
    )
    missing_remix_cellfields = np.array(set(remix_cellfields) - set(adata.obs.columns.values))
    print("MISSING CORPORA FIELDS:", missing_remix_cellfields)


# Process all
ad = sc.read_h5ad("huPr_Pd_all.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("huPr_Pd_all-curated.h5ad", compression="gzip")

rad = remix(ad, title="Normal Human Prostate Cell Atlas (All Cells)")
print_summary(rad)
rad.write("huPr_Pd_all-remixed.h5ad", compression="gzip")

# Process epithelial
ad = sc.read_h5ad("huPr_Pd_epi.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("huPr_Pd_epi-curated.h5ad", compression="gzip")

rad = remix(ad, title="Normal Human Prostate Cell Atlas (Epithelial Cells)")
print_summary(rad)
rad.write("huPr_Pd_epi-remixed.h5ad", compression="gzip")

# Process fibromuscular stromal cells
ad = sc.read_h5ad("huPr_Pd_fmst.h5ad")
basic_curation(ad)
print_summary(ad)
ad.write("huPr_Pd_fmst-curated.h5ad", compression="gzip")

rad = remix(ad, title="Normal Human Prostate Cell Atlas (Fibromuscular Stromal Cells)")
print_summary(rad)
rad.write("huPr_Pd_fmst-remixed.h5ad", compression="gzip")
