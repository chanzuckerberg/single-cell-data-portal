/* eslint-disable @typescript-eslint/camelcase */
/* eslint-disable sonarjs/no-duplicate-string */
import { DATASET_ASSET_FORMAT, DATASET_ASSET_TYPE, Project } from "../entities";

// Dataset URL example:
// "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/"

export const PROJECTS: Project[] = [
  {
    datasets: [
      {
        dataset_assets: [
          {
            file_name: "krasnow_lab_human_lung_cell_atlas_10x-1-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id: "krasnow_lab_human_lung_cell_atlas_10x-1-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name: "krasnow_lab_human_lung_cell_atlas_10x-1-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id: "krasnow_lab_human_lung_cell_atlas_10x-1-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id: "Krasnow Lab Human Lung Cell Atlas, 10X",
        title: "Krasnow Lab Human Lung Cell Atlas, 10X",
      },
      {
        dataset_assets: [
          {
            file_name: "krasnow_lab_human_lung_cell_atlas_smartseq2-2-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id: "krasnow_lab_human_lung_cell_atlas_smartseq2-2-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name: "krasnow_lab_human_lung_cell_atlas_smartseq2-2-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id: "krasnow_lab_human_lung_cell_atlas_smartseq2-2-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id: "Krasnow Lab Human Lung Cell Atlas, Smart-seq2",
        title: "Krasnow Lab Human Lung Cell Atlas, Smart-seq2",
      },
    ],
    links: [
      {
        name: "Krasnow Lab",
        url: "http://cmgm-new.stanford.edu/krasnow/",
      },
      {
        name: "HLCA website",
        url: "https://github.com/krasnowlab/hlca",
      },
    ],
  },
  {
    datasets: [
      {
        dataset_assets: [
          {
            file_name: "human_cell_landscape-3-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id: "human_cell_landscape-3-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name: "human_cell_landscape-3-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id: "human_cell_landscape-3-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id: "Human Cell Landscape",
        title: "Human Cell Landscape",
      },
    ],
    links: [
      {
        name: "Guo Lab",
        url: "https://person.zju.edu.cn/en/ggj",
      },
      {
        name: "HCL website",
        url: "http://bis.zju.edu.cn/HCL/",
      },
    ],
  },
  {
    datasets: [
      {
        dataset_assets: [
          {
            file_name:
              "human_fetal_liver_single_cell_transcriptome-13-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id: "human_fetal_liver_single_cell_transcriptome-13-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name: "human_fetal_liver_single_cell_transcriptome-13-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id: "human_fetal_liver_single_cell_transcriptome-13-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id: "Human fetal liver single cell transcriptome data",
        title: "Human fetal liver single cell transcriptome data",
      },
      {
        dataset_assets: [
          {
            file_name: "cell_atlas_of_thymic_development-14-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id: "cell_atlas_of_thymic_development-14-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name: "cell_atlas_of_thymic_development-14-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id: "cell_atlas_of_thymic_development-14-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id:
          "A cell atlas of human thymic development defines T cell repertoire formation",
        title:
          "A cell atlas of human thymic development defines T cell repertoire formation",
      },
    ],
    links: [
      {
        name: "E-MTAB-7407",
        url: "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7407/",
      },
      {
        name: "covid19cellatlas.org",
        url: "https://www.covid19cellatlas.org/",
      },
    ],
  },
  {
    datasets: [
      {
        dataset_assets: [
          {
            file_name:
              "cellular_census_of_human_lungs_alveoli_and_parenchyma-15-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id:
              "cellular_census_of_human_lungs_alveoli_and_parenchyma-15-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name:
              "cellular_census_of_human_lungs_alveoli_and_parenchyma-15-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id:
              "cellular_census_of_human_lungs_alveoli_and_parenchyma-15-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - parenchyma",
        title:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - parenchyma",
      },
      {
        dataset_assets: [
          {
            file_name: "cellular_census_of_human_lungs_nasal-16-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id: "cellular_census_of_human_lungs_nasal-16-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name: "cellular_census_of_human_lungs_nasal-16-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id: "cellular_census_of_human_lungs_nasal-16-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - nasal",
        title:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - nasal",
      },
      {
        dataset_assets: [
          {
            file_name: "cellular_census_of_human_lungs_bronchi-17-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id: "cellular_census_of_human_lungs_bronchi-17-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name: "cellular_census_of_human_lungs_bronchi-17-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id: "cellular_census_of_human_lungs_bronchi-17-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - bronchi",
        title:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - bronchi",
      },
      {
        dataset_assets: [
          {
            file_name:
              "ischaemic_sensitivity_of_human_tissue_by_single_cell_RNA_seq_lung-18-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id:
              "ischaemic_sensitivity_of_human_tissue_by_single_cell_RNA_seq_lung-18-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name:
              "ischaemic_sensitivity_of_human_tissue_by_single_cell_RNA_seq_lung-18-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id:
              "ischaemic_sensitivity_of_human_tissue_by_single_cell_RNA_seq_lung-18-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id:
          "Ischaemic sensitivity of human tissue by single cell RNA seq - lung",
        title:
          "Ischaemic sensitivity of human tissue by single cell RNA seq - lung",
      },
    ],
    links: [
      {
        name: "asthma.cellgeni.sanger.ac.uk",
        url: "https://asthma.cellgeni.sanger.ac.uk/",
      },
      {
        name: "covid19cellatlas.org",
        url: "https://www.covid19cellatlas.org/",
      },
    ],
  },
  {
    datasets: [
      {
        dataset_assets: [
          {
            file_name:
              "ischaemic_sensitivity_of_human_tissue_by_single_cell_RNA_seq_spleen-19-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id:
              "ischaemic_sensitivity_of_human_tissue_by_single_cell_RNA_seq_spleen-19-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name:
              "ischaemic_sensitivity_of_human_tissue_by_single_cell_RNA_seq_spleen-19-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id:
              "ischaemic_sensitivity_of_human_tissue_by_single_cell_RNA_seq_spleen-19-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id:
          "Ischaemic sensitivity of human tissue by single cell RNA seq - spleen",
        title:
          "Ischaemic sensitivity of human tissue by single cell RNA seq - spleen",
      },
    ],
    links: [
      {
        name: "HCA",
        url:
          "https://data.humancellatlas.org/explore/projects/c4077b3c-5c98-4d26-a614-246d12c2e5d7",
      },
      {
        name: "covid19cellatlas.org",
        url: "https://www.covid19cellatlas.org/",
      },
    ],
  },
  {
    datasets: [
      {
        dataset_assets: [
          {
            file_name:
              "Allergic_inflammatory_memory_in_human_respiratory_epithelial_progenitor_cells_scraping-9-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id:
              "Allergic_inflammatory_memory_in_human_respiratory_epithelial_progenitor_cells_scraping-9-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name:
              "Allergic_inflammatory_memory_in_human_respiratory_epithelial_progenitor_cells_scraping-9-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id:
              "Allergic_inflammatory_memory_in_human_respiratory_epithelial_progenitor_cells_scraping-9-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id:
          "Allergic inflammatory memory in human respiratory epithelial progenitor cells - nasal scrapings",
        title:
          "Allergic inflammatory memory in human respiratory epithelial progenitor cells - nasal scrapings",
      },
      {
        dataset_assets: [
          {
            file_name:
              "Allergic_inflammatory_memory_in_human_respiratory_epithelial_progenitor_cells_epithelial-10-original",
            format: DATASET_ASSET_FORMAT.CXG,
            id:
              "Allergic_inflammatory_memory_in_human_respiratory_epithelial_progenitor_cells_epithelial-10-original",
            s3_uri: "abc-original",
            type: DATASET_ASSET_TYPE.ORIGINAL,
          },
          {
            file_name:
              "Allergic_inflammatory_memory_in_human_respiratory_epithelial_progenitor_cells_epithelial-10-remix",
            format: DATASET_ASSET_FORMAT.H5AD,
            id:
              "Allergic_inflammatory_memory_in_human_respiratory_epithelial_progenitor_cells_epithelial-10-remix",
            s3_uri: "abc-remix",
            type: DATASET_ASSET_TYPE.REMIX,
          },
        ],
        id:
          "Allergic inflammatory memory in human respiratory epithelial progenitor cells - epithelial cells",
        title:
          "Allergic inflammatory memory in human respiratory epithelial progenitor cells - epithelial cells",
      },
    ],
    links: [
      {
        name: "Single Cell Portal",
        url:
          "https://singlecell.broadinstitute.org/single_cell/study/SCP253/allergic-inflammatory-memory-in-human-respiratory-epithelial-progenitor-cells?scpbr=the-alexandria-project",
      },
    ],
  },
];
