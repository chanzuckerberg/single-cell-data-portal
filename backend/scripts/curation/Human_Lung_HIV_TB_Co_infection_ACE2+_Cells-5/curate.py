import anndata
import pandas as pd

import utils.hgnc
import utils.ontology


def basic_curation(adata):

    adata.uns["organism_ontology"] = "NCBITaxon:9606"
    adata.uns["organism"] = utils.ontology.get_ontology_label(adata.uns["organism_ontology"])

    adata.uns["title"] = "Human Lung HIV-TB Co-infection ACE2+ Cells"
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

    disease_ontology_map = {
        "MONDO_0018076|MONDO_0005109": "MONDO:0018076|MONDO:0005109",
        "MONDO_0018076": "MONDO:0018076",
        "MONDO_0004992": "MONDO:0004992",
    }
    disease_map = {
        "MONDO_0018076|MONDO_0005109": "tuberculosis|HIV infectious disease",
        "MONDO_0018076": "tuberculosis",
        "MONDO_0004992": "cancer",
    }
    adata.obs["disease_ontology"] = adata.obs["disease"].replace(disease_ontology_map, inplace=False)
    adata.obs["disease"].replace(disease_map, inplace=True)

    ct_ontology_map = {}
    for cl_term in adata.obs["cell_type"].unique():
        cl_term_as_curie = cl_term.replace("_", ":")
        ct_ontology_map[cl_term] = cl_term_as_curie
    ct_map = {k: utils.ontology.get_ontology_label(v) for k, v in ct_ontology_map.items()}
    adata.obs["cell_type_ontology"] = adata.obs["cell_type"].replace(ct_ontology_map, inplace=False)
    adata.obs["cell_type"] = adata.obs["cell_type"].replace(ct_map, inplace=False)

    adata.obs["ethnicity"] = ""
    adata.obs["ethnicity_ontology"] = ""
    adata.obs["development_stage_ontology"] = "HsapDv:0000087"
    adata.obs["development_stage"] = utils.ontology.get_ontology_label(adata.obs["development_stage_ontology"][0])

    adata.uns["layer_descriptions"] = {"X": "raw"}

    upgraded_var_index = utils.hgnc.get_upgraded_var_index(adata.var)
    merged_df = pd.DataFrame(adata.X, index=adata.obs.index, columns=upgraded_var_index).sum(
        axis=1, level=0, skipna=False
    )

    remix_adata = anndata.AnnData(
        X=merged_df.to_numpy(),
        obs=adata.obs,
        var=merged_df.columns.to_frame(name="hgnc_gene_symbol"),
        uns=adata.uns,
        obsm=adata.obsm,
        varm=adata.varm,
    )

    return remix_adata


def create_adata():
    X_epi = pd.read_csv("Human_lung_epithelial_cell_raw_counts.txt.gz", sep="\t", header=0, index_col=0).T
    X_nonepi = pd.read_csv("Human_lung_nonepithelial_ACE2_TMPRSS2_raw_counts.txt.gz", sep="\t", header=0, index_col=0).T
    X = pd.concat([X_epi, X_nonepi])

    umap = pd.read_csv("AHRI_lung_umap.csv", sep=",", header=0, skiprows=[1], index_col=0)
    umap = umap.reindex(X.index)

    metadata = pd.read_csv("alexandria_structured_metadata.txt", sep="\t", header=0, index_col=0, skiprows=[1])
    metadata = metadata.reindex(X.index)

    # Remove direct personal identifiers
    # The geography refers to a specific city. DOC and DOB are full dates.
    metadata = metadata.drop(columns=["DOC", "DOB", "geographical_region", "geographical_region__ontology_label"])

    # Per submitter "ACE2_TMPRSS2 should be 'negative' where it is missing"
    metadata["ACE2_TMPRSS2_double_positive"].fillna("Negative", inplace=True)
    metadata = metadata.drop(columns=["CellID"])

    adata = anndata.AnnData(X=X, obs=metadata, obsm={"X_umap": umap.to_numpy()})
    adata.obs["donor_id"] = adata.obs["donor_id"].astype("category")
    adata.obs["biosample_id"] = adata.obs["biosample_id"].astype("category")

    return adata


def main():
    curated_filename = "Human_Lung_HIV_TB_Co_infection_ACE2+_Cells-5-curated.h5ad"
    remixed_filename = "Human_Lung_HIV_TB_Co_infection_ACE2+_Cells-5-remixed.h5ad"

    adata = create_adata()
    basic_curation(adata)
    adata.write(curated_filename, compression="gzip")
    remix_adata = remix(adata)
    remix_adata.write(remixed_filename, compression="gzip")


main()
