import anndata
import numpy as np
import pandas as pd
import scanpy.api as sc

import utils.hgnc
import utils.ontology


def basic_curation(adata):
    adata.obs["assay"] = "Microwell-seq"
    adata.obs["assay_ontology"] = ""
    adata.obs["disease_ontology"] = "PATO:0000461"
    adata.obs["disease"] = utils.ontology.get_ontology_label("PATO:0000461")

    adata.uns["organism_ontology"] = "NCBITaxon:9606"
    adata.uns["organism"] = utils.ontology.get_ontology_label("NCBITaxon:9606")
    adata.uns["title"] = "Construction of a human cell landscape at single-cell level"
    adata.uns["contributors"] = [
        {"name": "Huiyu Sun", "institution": "Zhejiang University School of Medicine", "email": "sunhuiyu@zju.edu.cn",},
        {"name": "Guoji Guo", "institution": "Zhejiang University School of Medicine", "email": "ggj@zju.edu.cn",},
    ]

    adata.uns["publication_doi"] = "https://doi.org/10.1038/s41586-020-2157-4"

    adata.uns["project_name"] = adata.uns["title"]
    adata.uns["project_description"] = (
        "Single-cell analysis is a valuable tool for dissecting cellular heterogeneity in complex systems. "
        "However, a comprehensive single-cell atlas has not been achieved for humans. Here we use single-cell "
        "mRNA sequencing to determine the cell-type composition of all major human organs and construct a scheme "
        "for the human cell landscape (HCL). We have uncovered a single-cell hierarchy for many tissues that have "
        "not been well characterized. We established a 'single-cell HCL analysis' pipeline that helps to define "
        "human cell identity. Finally, we performed a single-cell comparative analysis of landscapes from human "
        "and mouse to identify conserved genetic networks. We found that stem and progenitor cells exhibit "
        "strong transcriptomic stochasticity, whereas diferentiated cells are more distinct. Our results provide a"
        "useful resource for the study of human biology."
    )
    adata.uns["project_protocol_links"] = []
    adata.uns["project_raw_data_links"] = ["https://figshare.com/articles/HCL_DGE_Data/7235471"]
    adata.uns["project_other_links"] = [
        "https://db.cngb.org/HCL/",
        "https://github.com/ggjlab/HCL/",
    ]


def remix(adata):

    # Handle tissue. This one has lots of them, and the original name collides with the corpora name
    adata.obs.rename(columns={"tissue": "original_tissue"}, inplace=True)
    tissue_ontology_map = {
        "AdultLung": "UBERON:0002048",
        "FetalIntestine": "UBERON:0000160",
        "AdultAdrenalGland": "UBERON:0002369",
        "AdultKidney": "UBERON:0002113",
        "FetalKidney": "UBERON:0002113",
        "AdultPleura": "UBERON:0000977",
        "FetalPancreas": "UBERON:0001264",
        "FetalMuscle": "UBERON:0001630",
        "FetalLiver": "UBERON:0002107",
        "AdultPeripheralBlood": "UBERON:0000178",
        "AdultTransverseColon": "UBERON:0001157",
        "CordBloodCD34P": "UBERON:0012168",
        "AdultSpleen": "UBERON:0002106",
        "AdultStomach": "UBERON:0000945",
        "FetalAdrenalGland": "UBERON:0002369",
        "FetalBrain": "UBERON:0000955",
        "FetalMaleGonad": "UBERON:0000473",
        "AdultOmentum": "UBERON:0003688",
        "AdultThyroid": "UBERON:0002046",
        "AdultEsophagus": "UBERON:0001043",
        "AdultLiver": "UBERON:0002107",
        "AdultTrachea": "UBERON:0003126",
        "ChorionicVillus": "UBERON:0007106",
        "AdultGallbladder": "UBERON:0002110",
        "AdultPancreas": "UBERON:0001264",
        "AdultArtery": "UBERON:0001637",
        "FetalLung": "UBERON:0002048",
        "Placenta": "UBERON:0001987",
        "AdultTemporalLobe": "UBERON:0001871",
        "AdultBladder": "UBERON:0018707",
        "AdultBoneMarrow": "UBERON:0002371",
        "AdultCervix": "UBERON:0000002",
        "FetalHeart": "UBERON:0000948",
        "FetalStomach": "UBERON:0000945",
        "AdultMuscle": "UBERON:0001630",
        "AdultUterus": "UBERON:0000995",
        "AdultCerebellum": "UBERON:0002037",
        "FetalSkin": "UBERON:0002097",
        "FetalFemaleGonad": "UBERON:0000992",
        "CordBlood": "UBERON:0012168",
        "AdultFallopiantube": "UBERON:0003889",
        "FetalRib": "UBERON:0002228",
        "FetalSpinalCord": "UBERON:0002240",
        "NeonatalAdrenalGland": "UBERON:0002369",
        "AdultRectum": "UBERON:0001052",
        "AdultJeJunum": "UBERON:0002115",
        "FetalCalvaria": "UBERON:0004339",
        "AdultDuodenum": "UBERON:0002114",
        "FetalThymus": "UBERON:0002370",
        "AdultEpityphlon": "UBERON:0001154",
        "AdultIleum": "UBERON:0002116",
        "AdultSigmoidColon": "UBERON:0001159",
        "AdultHeart": "UBERON:0000948",
        "AdultProstate": "UBERON:0002367",
        "AdultUreter": "UBERON:0000056",
        "AdultAscendingColon": "UBERON:0001156",
        "FetalEyes": "UBERON:0000970",
        "HESC": "",
        "AdultAdipose": "UBERON:0001013",
    }
    tissue_map = {k: utils.ontology.get_ontology_label(v) for k, v in tissue_ontology_map.items()}
    tissue_map["HESC"] = "HESC"
    adata.obs["tissue_ontology"] = adata.obs["original_tissue"].replace(tissue_ontology_map, inplace=False)
    adata.obs["tissue"] = adata.obs["original_tissue"].replace(tissue_map, inplace=False)

    adata.obs["cell_type_ontology"] = ""
    adata.obs["cell_type"] = ""
    adata.obs["ethnicity_ontology"] = "HANCESTRO:0027"
    adata.obs["ethnicity"] = utils.ontology.get_ontology_label("HANCESTRO:0027")
    adata.obs["sex"] = "unknown"

    development_stage_ontology_map = {}
    for k in tissue_ontology_map:
        if k.startswith("Adult") or k in ("CordBlood", "Placenta", "ChorionicVillus", "CordBloodCD34P",):
            development_stage_ontology_map[k] = "HsapDv:0000087"
        elif k.startswith("Fetal"):
            development_stage_ontology_map[k] = "HsapDv:0000037"
        elif k.startswith("Neonatal"):
            development_stage_ontology_map[k] = "HsapDv:0000082"
        elif k == "HESC":
            development_stage_ontology_map[k] = "HsapDv:0000002"
    development_stage_map = {k: utils.ontology.get_ontology_label(v) for k, v in development_stage_ontology_map.items()}
    adata.obs["development_stage_ontology"] = adata.obs["original_tissue"].replace(
        development_stage_ontology_map, inplace=False
    )
    adata.obs["development_stage"] = adata.obs["original_tissue"].replace(development_stage_map, inplace=False)

    adata.uns["layer_descriptions"] = {"X": "log1p CPM"}

    upgraded_var_index = utils.hgnc.get_upgraded_var_index(adata.var)
    merged_df = pd.DataFrame(np.expm1(adata.X), index=adata.obs.index, columns=upgraded_var_index).sum(
        axis=1, level=0, skipna=False
    )

    remix_adata = anndata.AnnData(
        X=np.log1p(merged_df.to_numpy()),
        obs=adata.obs,
        var=merged_df.columns.to_frame(name="hgnc_gene_symbol"),
        uns=adata.uns,
        obsm=adata.obsm,
        varm=adata.varm,
    )

    return remix_adata


def main():
    original_filename = "human_cell_landscape-3-original.h5ad"
    curated_filename = "human_cell_landscape-3-curated.h5ad"
    remixed_filename = "human_cell_landscape-3-remixed.h5ad"

    # Read raw, X has most of the genes filtered out
    adata = sc.read_h5ad(original_filename).raw.to_adata()
    basic_curation(adata)
    adata.write(curated_filename, compression="gzip")
    remix_adata = remix(adata)
    remix_adata.write(remixed_filename, compression="gzip")


main()
