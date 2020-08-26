import anndata
import pandas as pd

import utils.ontology


def basic_curation(adata):

    adata.uns["organism_ontology"] = "NCBITaxon:9541"
    adata.uns["organism"] = utils.ontology.get_ontology_label(adata.uns["organism_ontology"])
    adata.uns["title"] = "Epithelial Cells in NHP mTB Granuloma and Uninvolved Lung"
    adata.uns["contributors"] = [
        {
            "name": "Jose Ordovas-Montanes",
            "institution": "Broad Institute of MIT and Harvard",
            "email": "jose.ordovas-montanes@childrens.harvard.edu",
        },
        {"name": "Sarah Nyquist", "institution": " Massachusetts Institute of Technology", "email": "nyquist@mit.edu"},
    ]

    adata.uns["publication_doi"] = "https://doi.org/10.1016/j.cell.2020.04.035"

    adata.uns["project_name"] = adata.uns["title"]
    adata.uns["project_description"] = (
        "There is pressing urgency to understand the pathogenesis of the severe acute respiratory syndrome "
        "coronavirus clade 2 (SARS-CoV-2), which causes the disease COVID-19. SARS-CoV-2 spike (S) protein "
        "binds angiotensin-converting enzyme 2 (ACE2), and in concert with host proteases, principally "
        "transmembrane serine protease 2 (TMPRSS2), promotes cellular entry. The cell subsets targeted by "
        "SARS-CoV-2 in host tissues and the factors that regulate ACE2 expression remain unknown. Here, we "
        "leverage human, non-human primate, and mouse single-cell RNA-sequencing (scRNA-seq) datasets "
        "across health and disease to uncover putative targets of SARS-CoV-2 among tissue-resident cell "
        "subsets. We identify ACE2 and TMPRSS2 co-expressing cells within lung type II pneumocytes, ileal "
        "absorptive enterocytes, and nasal goblet secretory cells. Strikingly, we discovered that ACE2 is "
        "a human interferon-stimulated gene (ISG) in vitro using airway epithelial cells and extend our "
        "findings to in vivo viral infections. Our data suggest that SARS-CoV-2 could exploit "
        "species-specific interferon-driven upregulation of ACE2, a tissue-protective mediator during lung "
        "injury, to enhance infection."
    )
    adata.uns["project_protocol_links"] = []
    adata.uns["project_raw_data_links"] = ["http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148829"]
    adata.uns["project_other_links"] = [
        "https://singlecell.broadinstitute.org/single_cell?scpbr=the-alexandria-project",
        "https://singlecell.broadinstitute.org/single_cell/covid19",
    ]


def remix(adata):
    adata.obs["assay_ontology"] = "EFO:0008919"
    adata.obs["assay"] = utils.ontology.get_ontology_label(adata.obs["assay_ontology"][0])

    tissue_ontology_map = {}
    for uberon in adata.obs["organ"].unique():
        uberon_as_curie = uberon.replace("_", ":")
        tissue_ontology_map[uberon] = uberon_as_curie
    tissue_map = {k: utils.ontology.get_ontology_label(v) for k, v in tissue_ontology_map.items()}
    adata.obs["tissue_ontology"] = adata.obs["organ"].replace(tissue_ontology_map, inplace=False)
    adata.obs["tissue"] = adata.obs["organ"].replace(tissue_map, inplace=False)

    adata.obs["disease_ontology"] = "MONDO:0018076"
    adata.obs["disease"] = utils.ontology.get_ontology_label(adata.obs["disease_ontology"][0])

    ct_ontology_map = {}
    for cl_term in adata.obs["cell_type"].unique():
        cl_term_as_curie = cl_term.replace("_", ":")
        ct_ontology_map[cl_term] = cl_term_as_curie
    ct_map = {k: utils.ontology.get_ontology_label(v) for k, v in ct_ontology_map.items()}
    adata.obs["cell_type_ontology"] = adata.obs["cell_type"].replace(ct_ontology_map, inplace=False)
    adata.obs["cell_type"] = adata.obs["cell_type"].replace(ct_map, inplace=False)

    adata.obs["ethnicity"] = "na"
    adata.obs["ethnicity_ontology"] = ""
    adata.obs["development_stage_ontology"] = "EFO:0001272"
    adata.obs["development_stage"] = utils.ontology.get_ontology_label(adata.obs["development_stage_ontology"][0])

    adata.uns["layer_descriptions"] = {
        "X": "seurat lognormalize, scale.factor = 10000",
        "raw.X": "raw",
    }


def create_adata():

    X = pd.read_csv("epithelial_cells_for_scp.txt", sep="\t", header=0, index_col=0).T
    raw_X = pd.read_csv("all_epi_counts.txt", sep="\t", header=0, index_col=0)
    assert X.shape == raw_X.shape
    assert (X.columns == raw_X.columns).all()
    assert (X.index == raw_X.index).all()

    umap = pd.read_csv("epi_cells_umap.txt", sep="\t", header=0, skiprows=[1], index_col=0)
    umap = umap.reindex(X.index)
    metadata = pd.read_csv(
        "epithelial_metadata_alexandria.txt",
        sep="\t",
        header=0,
        index_col=0,
        skiprows=[1],
    )
    metadata = metadata.reindex(X.index)

    metadata = metadata.drop(columns=["CellID"])  # This is already the index

    adata = anndata.AnnData(X=X, obs=metadata, obsm={"X_umap": umap.to_numpy()})
    adata.raw = anndata.AnnData(X=raw_X)

    adata.obs["donor_id"] = adata.obs["donor_id"].astype("category")

    return adata


def main():
    curated_filename = "Epithelial_Cells_in_NHP_mTB_Granuloma_and_Uninvolved_Lung-4-curated.h5ad"
    remixed_filename = "Epithelial_Cells_in_NHP_mTB_Granuloma_and_Uninvolved_Lung-4-remixed.h5ad"
    adata = create_adata()
    basic_curation(adata)
    adata.write(curated_filename, compression="gzip")
    remix(adata)
    adata.write(remixed_filename, compression="gzip")


main()
