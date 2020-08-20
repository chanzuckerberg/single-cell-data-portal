/* eslint-disable @typescript-eslint/camelcase */
/* eslint-disable sonarjs/no-duplicate-string */
import { DATASET_ASSET_TYPE, Project } from "../entities";

// Dataset URL example:
// "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/"

export const PROJECTS: Project[] = [
  {
    datasets: [
      {
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id: "krasnow_lab_human_lung_cell_atlas_10x-1-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id: "Krasnow Lab Human Lung Cell Atlas, 10X",
        name: "Krasnow Lab Human Lung Cell Atlas, 10X",
      },
      {
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id: "krasnow_lab_human_lung_cell_atlas_smartseq2-2-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id: "Krasnow Lab Human Lung Cell Atlas, Smart-seq2",
        name: "Krasnow Lab Human Lung Cell Atlas, Smart-seq2",
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
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id: "human_cell_landscape-3-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id: "Human Cell Landscape",
        name: "Human Cell Landscape",
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
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id: "human_fetal_liver_single_cell_transcriptome-13-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id: "Human fetal liver single cell transcriptome data",
        name: "Human fetal liver single cell transcriptome data",
      },
      {
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id: "cell_atlas_of_thymic_development-14-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id:
          "A cell atlas of human thymic development defines T cell repertoire formation",
        name:
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
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id:
              "cellular_census_of_human_lungs_alveoli_and_parenchyma-15-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - parenchyma",
        name:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - parenchyma",
      },
      {
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id: "cellular_census_of_human_lungs_nasal-16-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - nasal",
        name:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - nasal",
      },
      {
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id: "cellular_census_of_human_lungs_bronchi-17-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - bronchi",
        name:
          "A cellular census of human lungs identifies novel cell states in health and in asthma - bronchi",
      },
      {
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id:
              "ischaemic_sensitivity_of_human_tissue_by_single_cell_RNA_seq_lung-18-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id:
          "Ischaemic sensitivity of human tissue by single cell RNA seq - lung",
        name:
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
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id:
              "ischaemic_sensitivity_of_human_tissue_by_single_cell_RNA_seq_spleen-19-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id:
          "Ischaemic sensitivity of human tissue by single cell RNA seq - spleen",
        name:
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
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id:
              "Allergic_inflammatory_memory_in_human_respiratory_epithelial_progenitor_cells_scraping-9-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id:
          "Allergic inflammatory memory in human respiratory epithelial progenitor cells - nasal scrapings",
        name:
          "Allergic inflammatory memory in human respiratory epithelial progenitor cells - nasal scrapings",
      },
      {
        dataset_deployments: [
          {
            environment: DATASET_ASSET_TYPE.REMIX,
            id:
              "Allergic_inflammatory_memory_in_human_respiratory_epithelial_progenitor_cells_epithelial-10-remix",
            link:
              "https://cellxgene.cziscience.com/d/krasnow_lab_human_lung_cell_atlas_10x-1.cxg/",
          },
        ],
        id:
          "Allergic inflammatory memory in human respiratory epithelial progenitor cells - epithelial cells",
        name:
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
