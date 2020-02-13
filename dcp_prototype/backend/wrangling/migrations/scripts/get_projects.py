from dcplib.etl import DSSExtractor

if __name__ == "__main__":
    projects = [
        "HumanMousePancreas",
        "1M Neurons",
        "MouseGastrulationAtlas",
        "HumanTissueTcellActivation",
        "WongAdultRetina",
        "HPSI human cerebral organoids",
        "BM_PC",
        "Multiplexed scRNA-seq with barcoded antibodies",
        "1M Immune Cells",
        "Drop-seq, DroNc-seq, Fluidigm C1 Comparison",
        "SingleCellLiverLandscape",
        "HDCA-Sweden-10x",
        "TissueStability",
        "Mouse Melanoma",
        "CD4+ cytotoxic T lymphocytes",
        "Human Hematopoietic Profiling",
        "Fetal/Maternal Interface",
        "Kidney biopsy scRNA-seq",
        "Reprogrammed_Dendritic_Cells",
        "Single cell transcriptome analysis of human pancreas",
        "Single cell RNAseq characterization of cell types produced over time in an in vitro model of human inhibitory interneuron differentiation",
        "Healthy and type 2 diabetes pancreas",
        "KidneySingleCellAtlas",
        "HumanColonicMesenchymeIBD",
        "scRNAseqSystemicComparison",
        "Tabula Muris",
        "Diabetic Nephropathy snRNA-seq",
        "Mouse Endoderm Project",
        "snRNA-seq_for_human_retina"]
    for project in projects:
        query = {"query": {
            "bool": {"must": [{"term": {"files.project_json.project_core.project_short_name": f"{project}"}}]}}}
        project = project.replace(" ", "-")
        print(project)
        DSSExtractor(staging_directory=f"./hold_data/{project}").extract(query=query, max_workers=50)
