import anndata
import pandas as pd

import utils.hgnc
import utils.ontology


def basic_curation(adata):

    adata.uns["organism_ontology"] = "NCBITaxon:9606"
    adata.uns["organism"] = utils.ontology.get_ontology_label(adata.uns["organism_ontology"])

    adata.uns["title"] = "ACE2 and TMPRSS2 expression in human non-inflamed terminal ileum"
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
    adata.obs["assay_ontology"] = "EFO:0009899"
    adata.obs["assay"] = utils.ontology.get_ontology_label(adata.obs["assay_ontology"][0])
    adata.obs["tissue_ontology"] = "UBERON:0000331"
    adata.obs["tissue"] = utils.ontology.get_ontology_label(adata.obs["tissue_ontology"][0])
    adata.obs["disease_ontology"] = "MONDO:0001318"
    adata.obs["disease"] = utils.ontology.get_ontology_label(adata.obs["disease_ontology"][0])

    ct_ontology_map = {}
    for cl_term in adata.obs["cell_type"].unique():
        cl_term_as_curie = cl_term.replace("_", ":")
        ct_ontology_map[cl_term] = cl_term_as_curie
    ct_map = {k: utils.ontology.get_ontology_label(v) for k, v in ct_ontology_map.items()}
    adata.obs["cell_type_ontology"] = adata.obs["cell_type"].replace(ct_ontology_map, inplace=False)
    adata.obs["cell_type"] = adata.obs["cell_type"].replace(ct_map, inplace=False)

    adata.obs["ethnicity"] = ""
    adata.obs["ethnicity_ontology"] = ""
    adata.obs["development_stage_ontology"] = "HsapDv:0000081"
    adata.obs["development_stage"] = utils.ontology.get_ontology_label(adata.obs["development_stage_ontology"][0])

    adata.uns["layer_descriptions"] = {"raw.X": "raw", "X": "seurat lognormalize"}

    # Cannot fix the gene names because not all the expression values are available


def create_adatas():
    # The expression values are spread across two csvs, and one of them only has ACE2 and TMPRSS2
    X_1 = pd.read_csv(
        "SCPdata_seurat_normalized_expression_counts_subset_ACE2_TMPRSS2_(1).csv", sep=",", header=0, index_col=0
    ).T
    X_2 = pd.read_csv(
        "SCPdata_seurat_normalized_expression_counts_subset_absorptiveAndCryptentero.csv",
        sep=",",
        header=0,
        index_col=0,
    ).T
    X = pd.concat([X_1, X_2], copy=False)

    # Same with the raw counts
    rawX_1 = pd.read_csv("seurat_counts_subset_notEC_ACE2_TMPRSS2.csv", sep=",", header=0, index_col=0).T
    rawX_2 = pd.read_csv("seurat_counts_subset_absorptiveAndCryptentero.csv", sep=",", header=0, index_col=0).T
    rawX = pd.concat([rawX_1, rawX_2], copy=False)

    assert (X.index == rawX.index).all()
    assert (X.columns == rawX.columns).all()
    assert X.shape == rawX.shape

    # There are two tsnes, one for everything and one for just the epithelial cells
    tsne_all = pd.read_csv(
        "SCPdata_20200311_tsne_all_noninflamed.csv", sep=",", header=0, skiprows=[1], index_col=0, usecols=[0, 1, 2]
    )
    tsne_epith = pd.read_csv(
        "SCPdata_20200311_tsne_epith_noninflamed.csv", sep=",", header=0, skiprows=[1], index_col=0, usecols=[0, 1, 2]
    )

    metadata = pd.read_csv("SCPdata_20200311_prelim_metadata.csv", sep=",", header=0, index_col=0, skiprows=[1])

    # Subset by the tsne indices
    X_all = X.reindex(tsne_all.index, copy=True)
    rawX_all = rawX.reindex(tsne_all.index, copy=True)
    X_epith = X.reindex(tsne_epith.index, copy=True)
    rawX_epith = rawX.reindex(tsne_epith.index, copy=True)

    metadata_all = metadata.reindex(X_all.index, copy=True)
    metadata_epith = metadata.reindex(X_epith.index, copy=True)

    adata_all = anndata.AnnData(
        X=X_all, obs=metadata_all, obsm={"X_tsne": tsne_all.to_numpy()}, raw=anndata.AnnData(X=rawX_all)
    )
    adata_epith = anndata.AnnData(
        X=X_epith, obs=metadata_epith, obsm={"X_tsne": tsne_epith.to_numpy()}, raw=anndata.AnnData(X=rawX_epith)
    )

    return adata_all, adata_epith


def main():

    all_curated_fn = "ACE2_and_TMPRSS2_expression_in_human_non_inflamed_terminal_ileum-6-curated.h5ad"
    all_remixed_fn = "ACE2_and_TMPRSS2_expression_in_human_non_inflamed_terminal_ileum-6-remixed.h5ad"
    epi_curated_fn = "ACE2_and_TMPRSS2_expression_in_human_non_inflamed_terminal_ileum_epithelial-7-curated.h5ad"
    epi_remixed_fn = "ACE2_and_TMPRSS2_expression_in_human_non_inflamed_terminal_ileum_epithelial-7-remixed.h5ad"

    adata_all, adata_epi = create_adatas()

    basic_curation(adata_all)
    adata_all.write(all_curated_fn, compression="gzip")
    basic_curation(adata_epi)
    adata_epi.write(epi_curated_fn, compression="gzip")

    remix(adata_all)
    adata_all.write(all_remixed_fn, compression="gzip")
    remix(adata_epi)
    adata_epi.write(epi_remixed_fn, compression="gzip")


main()
