/* eslint-disable sonarjs/no-duplicate-string */

import { Dataset } from "src/common/entities";

const datasets = [
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 3961,
    collection_id: "b52eb423-5d0d-4645-b217-e1c6d38b2e72",
    collection_visibility: "PUBLIC",
    created_at: 1617746933.623842,
    dataset_assets: [
      {
        created_at: 1617747599.150022,
        dataset_id: "1009f384-b12d-448e-ba9f-1b7d2ecfbb4e",
        filename: "local.h5ad",
        filetype: "H5AD",
        id: "d925fb89-8b05-4aff-b9c3-00c810cb71bd",
        s3_uri:
          "s3://corpora-data-prod/1009f384-b12d-448e-ba9f-1b7d2ecfbb4e/local.h5ad",
        type: "REMIX",
        updated_at: 1617747599.150026,
        user_submitted: true,
      },
      {
        created_at: 1617747639.017621,
        dataset_id: "1009f384-b12d-448e-ba9f-1b7d2ecfbb4e",
        filename: "local.rds",
        filetype: "RDS",
        id: "4d0462e2-621e-4d80-94c6-9577b2215373",
        s3_uri:
          "s3://corpora-data-prod/1009f384-b12d-448e-ba9f-1b7d2ecfbb4e/local.rds",
        type: "REMIX",
        updated_at: 1617747639.017628,
        user_submitted: true,
      },
    ],
    dataset_deployments: [
      {
        dataset_id: "1009f384-b12d-448e-ba9f-1b7d2ecfbb4e",
        url: "https://cellxgene.cziscience.com/e/1009f384-b12d-448e-ba9f-1b7d2ecfbb4e.cxg/",
      },
    ],
    development_stage: [
      {
        label: "40-year-old human stage",
        ontology_term_id: "HsapDv:0000134",
      },
      {
        label: "45-year-old human stage",
        ontology_term_id: "HsapDv:0000139",
      },
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "55-year-old human stage",
        ontology_term_id: "HsapDv:0000149",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "65-year-old human stage",
        ontology_term_id: "HsapDv:0000159",
      },
      {
        label: "70-year-old human stage",
        ontology_term_id: "HsapDv:0000164",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "1009f384-b12d-448e-ba9f-1b7d2ecfbb4e",
    is_valid: false,
    linked_genesets: [],
    name: "Neuronal \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    processing_status: {
      h5ad_status: "CONVERTED",
      cxg_status: "CONVERTED",
      rds_status: "CONVERTED",
      created_at: 1617746933.625826,
      dataset_id: "1009f384-b12d-448e-ba9f-1b7d2ecfbb4e",
      id: "4a17ba9b-ea54-42dd-8eee-2e20b6160961",
      processing_status: "SUCCESS",
      updated_at: 1617747639.074367,
      upload_progress: 1.0,
      upload_status: "UPLOADED",
      validation_status: "VALID",
    },
    published: true,
    revision: 0,
    sex: ["female", "male"],
    tissue: [
      {
        label: "apex of heart",
        ontology_term_id: "UBERON:0002098",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "interventricular septum",
        ontology_term_id: "UBERON:0002094",
      },
      {
        label: "left cardiac atrium",
        ontology_term_id: "UBERON:0002079",
      },
      {
        label: "right cardiac atrium",
        ontology_term_id: "UBERON:0002078",
      },
    ],
    updated_at: 1617755780.245373,
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 195395,
    collection_id: "b52eb423-5d0d-4645-b217-e1c6d38b2e72",
    collection_visibility: "PUBLIC",
    created_at: 1617746940.676077,
    dataset_assets: [
      {
        created_at: 1617748198.665936,
        dataset_id: "572f3f3e-d3e4-4d13-8e2b-88215e508481",
        filename: "local.h5ad",
        filetype: "H5AD",
        id: "ca52fdae-f092-4fbf-afb9-58cb7486a0ca",
        s3_uri:
          "s3://corpora-data-prod/572f3f3e-d3e4-4d13-8e2b-88215e508481/local.h5ad",
        type: "REMIX",
        updated_at: 1617748198.66594,
        user_submitted: true,
      },
      {
        created_at: 1617749879.638822,
        dataset_id: "572f3f3e-d3e4-4d13-8e2b-88215e508481",
        filename: "local.rds",
        filetype: "RDS",
        id: "a3a4ea58-da04-4fcf-8cc8-db9f47b4d58a",
        s3_uri:
          "s3://corpora-data-prod/572f3f3e-d3e4-4d13-8e2b-88215e508481/local.rds",
        type: "REMIX",
        updated_at: 1617749879.638827,
        user_submitted: true,
      },
    ],
    dataset_deployments: [
      {
        dataset_id: "572f3f3e-d3e4-4d13-8e2b-88215e508481",
        url: "https://cellxgene.cziscience.com/e/572f3f3e-d3e4-4d13-8e2b-88215e508481.cxg/",
      },
    ],
    development_stage: [
      {
        label: "40-year-old human stage",
        ontology_term_id: "HsapDv:0000134",
      },
      {
        label: "45-year-old human stage",
        ontology_term_id: "HsapDv:0000139",
      },
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "55-year-old human stage",
        ontology_term_id: "HsapDv:0000149",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "65-year-old human stage",
        ontology_term_id: "HsapDv:0000159",
      },
      {
        label: "70-year-old human stage",
        ontology_term_id: "HsapDv:0000164",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "572f3f3e-d3e4-4d13-8e2b-88215e508481",
    is_valid: false,
    linked_genesets: [],
    name: "Vascular \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    processing_status: {
      h5ad_status: "CONVERTED",
      cxg_status: "CONVERTED",
      rds_status: "CONVERTED",
      created_at: 1617746940.678959,
      dataset_id: "572f3f3e-d3e4-4d13-8e2b-88215e508481",
      id: "1ec36b20-447a-4f5c-8f53-53e8ef4f0362",
      processing_status: "SUCCESS",
      updated_at: 1617749879.731402,
      upload_progress: 1.0,
      upload_status: "UPLOADED",
      validation_status: "VALID",
    },
    published: true,
    revision: 0,
    sex: ["female", "male"],
    tissue: [
      {
        label: "apex of heart",
        ontology_term_id: "UBERON:0002098",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "interventricular septum",
        ontology_term_id: "UBERON:0002094",
      },
      {
        label: "left cardiac atrium",
        ontology_term_id: "UBERON:0002079",
      },
      {
        label: "right cardiac atrium",
        ontology_term_id: "UBERON:0002078",
      },
    ],
    updated_at: 1617755780.245383,
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 125289,
    collection_id: "b52eb423-5d0d-4645-b217-e1c6d38b2e72",
    collection_visibility: "PUBLIC",
    created_at: 1617746880.071121,
    dataset_assets: [
      {
        created_at: 1617747955.808853,
        dataset_id: "78fd69d2-75e4-4207-819a-563139f273c6",
        filename: "local.h5ad",
        filetype: "H5AD",
        id: "dc201379-1e9e-4973-ad6f-8367e068f4d3",
        s3_uri:
          "s3://corpora-data-prod/78fd69d2-75e4-4207-819a-563139f273c6/local.h5ad",
        type: "REMIX",
        updated_at: 1617747955.808857,
        user_submitted: true,
      },
      {
        created_at: 1617749140.301681,
        dataset_id: "78fd69d2-75e4-4207-819a-563139f273c6",
        filename: "local.rds",
        filetype: "RDS",
        id: "178499a3-eaa0-4617-aee2-1e618cf17bda",
        s3_uri:
          "s3://corpora-data-prod/78fd69d2-75e4-4207-819a-563139f273c6/local.rds",
        type: "REMIX",
        updated_at: 1617749140.301688,
        user_submitted: true,
      },
    ],
    dataset_deployments: [
      {
        dataset_id: "78fd69d2-75e4-4207-819a-563139f273c6",
        url: "https://cellxgene.cziscience.com/e/78fd69d2-75e4-4207-819a-563139f273c6.cxg/",
      },
    ],
    development_stage: [
      {
        label: "40-year-old human stage",
        ontology_term_id: "HsapDv:0000134",
      },
      {
        label: "45-year-old human stage",
        ontology_term_id: "HsapDv:0000139",
      },
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "55-year-old human stage",
        ontology_term_id: "HsapDv:0000149",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "65-year-old human stage",
        ontology_term_id: "HsapDv:0000159",
      },
      {
        label: "70-year-old human stage",
        ontology_term_id: "HsapDv:0000164",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "78fd69d2-75e4-4207-819a-563139f273c6",
    is_valid: false,
    linked_genesets: [],
    name: "Ventricular cardiomyocytes \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    processing_status: {
      h5ad_status: "CONVERTED",
      cxg_status: "CONVERTED",
      rds_status: "CONVERTED",
      created_at: 1617746880.074091,
      dataset_id: "78fd69d2-75e4-4207-819a-563139f273c6",
      id: "846291c9-6256-43e0-9794-bbd218e2a8e6",
      processing_status: "SUCCESS",
      updated_at: 1617749140.359666,
      upload_progress: 1.0,
      upload_status: "UPLOADED",
      validation_status: "VALID",
    },
    published: true,
    revision: 0,
    sex: ["female", "male"],
    tissue: [
      {
        label: "apex of heart",
        ontology_term_id: "UBERON:0002098",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "interventricular septum",
        ontology_term_id: "UBERON:0002094",
      },
    ],
    updated_at: 1617755780.245388,
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 23483,
    collection_id: "b52eb423-5d0d-4645-b217-e1c6d38b2e72",
    collection_visibility: "PUBLIC",
    created_at: 1617746868.291982,
    dataset_assets: [
      {
        created_at: 1617747325.537067,
        dataset_id: "84f1a631-910b-4fbb-9f76-d915a07316d2",
        filename: "local.h5ad",
        filetype: "H5AD",
        id: "c5daa006-b246-42a1-ad2a-9fad66b414e8",
        s3_uri:
          "s3://corpora-data-prod/84f1a631-910b-4fbb-9f76-d915a07316d2/local.h5ad",
        type: "REMIX",
        updated_at: 1617747325.537073,
        user_submitted: true,
      },
      {
        created_at: 1617747541.115583,
        dataset_id: "84f1a631-910b-4fbb-9f76-d915a07316d2",
        filename: "local.rds",
        filetype: "RDS",
        id: "06440235-cff1-4242-b295-3fd361bf23f4",
        s3_uri:
          "s3://corpora-data-prod/84f1a631-910b-4fbb-9f76-d915a07316d2/local.rds",
        type: "REMIX",
        updated_at: 1617747541.115589,
        user_submitted: true,
      },
    ],
    dataset_deployments: [
      {
        dataset_id: "84f1a631-910b-4fbb-9f76-d915a07316d2",
        url: "https://cellxgene.cziscience.com/e/84f1a631-910b-4fbb-9f76-d915a07316d2.cxg/",
      },
    ],
    development_stage: [
      {
        label: "40-year-old human stage",
        ontology_term_id: "HsapDv:0000134",
      },
      {
        label: "45-year-old human stage",
        ontology_term_id: "HsapDv:0000139",
      },
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "55-year-old human stage",
        ontology_term_id: "HsapDv:0000149",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "65-year-old human stage",
        ontology_term_id: "HsapDv:0000159",
      },
      {
        label: "70-year-old human stage",
        ontology_term_id: "HsapDv:0000164",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "84f1a631-910b-4fbb-9f76-d915a07316d2",
    is_valid: false,
    linked_genesets: [],
    name: "Atrial cardiomyocytes \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    processing_status: {
      h5ad_status: "CONVERTED",
      cxg_status: "CONVERTED",
      rds_status: "CONVERTED",
      created_at: 1617746868.294075,
      dataset_id: "84f1a631-910b-4fbb-9f76-d915a07316d2",
      id: "a74da350-6d65-411c-9f56-685d3f67707a",
      processing_status: "SUCCESS",
      updated_at: 1617747541.203532,
      upload_progress: 1.0,
      upload_status: "UPLOADED",
      validation_status: "VALID",
    },
    published: true,
    revision: 0,
    sex: ["female", "male"],
    tissue: [
      {
        label: "left cardiac atrium",
        ontology_term_id: "UBERON:0002079",
      },
      {
        label: "right cardiac atrium",
        ontology_term_id: "UBERON:0002078",
      },
    ],
    updated_at: 1617755780.245391,
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 59341,
    collection_id: "b52eb423-5d0d-4645-b217-e1c6d38b2e72",
    collection_visibility: "PUBLIC",
    created_at: 1617746890.093482,
    dataset_assets: [
      {
        created_at: 1617747653.571758,
        dataset_id: "9d584fcb-a28a-4b91-a886-ceb66a88ef81",
        filename: "local.h5ad",
        filetype: "H5AD",
        id: "824875af-c989-43d9-9c07-d782fd06f5bc",
        s3_uri:
          "s3://corpora-data-prod/9d584fcb-a28a-4b91-a886-ceb66a88ef81/local.h5ad",
        type: "REMIX",
        updated_at: 1617747653.571765,
        user_submitted: true,
      },
      {
        created_at: 1617748145.260177,
        dataset_id: "9d584fcb-a28a-4b91-a886-ceb66a88ef81",
        filename: "local.rds",
        filetype: "RDS",
        id: "8ab73c42-8086-427f-b585-3ada3e7f55f3",
        s3_uri:
          "s3://corpora-data-prod/9d584fcb-a28a-4b91-a886-ceb66a88ef81/local.rds",
        type: "REMIX",
        updated_at: 1617748145.260182,
        user_submitted: true,
      },
    ],
    dataset_deployments: [
      {
        dataset_id: "9d584fcb-a28a-4b91-a886-ceb66a88ef81",
        updated_at: 1617747615.886036,
        url: "https://cellxgene.cziscience.com/e/9d584fcb-a28a-4b91-a886-ceb66a88ef81.cxg/",
      },
    ],
    development_stage: [
      {
        label: "40-year-old human stage",
        ontology_term_id: "HsapDv:0000134",
      },
      {
        label: "45-year-old human stage",
        ontology_term_id: "HsapDv:0000139",
      },
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "55-year-old human stage",
        ontology_term_id: "HsapDv:0000149",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "65-year-old human stage",
        ontology_term_id: "HsapDv:0000159",
      },
      {
        label: "70-year-old human stage",
        ontology_term_id: "HsapDv:0000164",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "9d584fcb-a28a-4b91-a886-ceb66a88ef81",
    is_valid: false,
    linked_genesets: [],
    name: "Fibroblasts \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    processing_status: {
      h5ad_status: "CONVERTED",
      cxg_status: "CONVERTED",
      rds_status: "CONVERTED",
      created_at: 1617746890.095562,
      dataset_id: "9d584fcb-a28a-4b91-a886-ceb66a88ef81",
      id: "73beb21b-33d3-4123-8bed-83ed5d830784",
      processing_status: "SUCCESS",
      updated_at: 1617748145.334961,
      upload_progress: 1.0,
      upload_status: "UPLOADED",
      validation_status: "VALID",
    },
    published: true,
    revision: 0,
    sex: ["female", "male"],
    tissue: [
      {
        label: "apex of heart",
        ontology_term_id: "UBERON:0002098",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "interventricular septum",
        ontology_term_id: "UBERON:0002094",
      },
      {
        label: "left cardiac atrium",
        ontology_term_id: "UBERON:0002079",
      },
      {
        label: "right cardiac atrium",
        ontology_term_id: "UBERON:0002078",
      },
    ],
    updated_at: 1617755780.245393,
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 486134,
    collection_id: "b52eb423-5d0d-4645-b217-e1c6d38b2e72",
    collection_visibility: "PUBLIC",
    created_at: 1617319340.173896,
    dataset_assets: [
      {
        created_at: 1617321217.25322,
        dataset_id: "d4e69e01-3ba2-4d6b-a15d-e7048f78f22e",
        filename: "local.h5ad",
        filetype: "H5AD",
        id: "dada2b31-236c-4585-95fc-05d91fba3bef",
        s3_uri:
          "s3://corpora-data-prod/d4e69e01-3ba2-4d6b-a15d-e7048f78f22e/local.h5ad",
        type: "REMIX",
        updated_at: 1617321217.253226,
        user_submitted: true,
      },
      {
        created_at: 1617323095.801804,
        dataset_id: "d4e69e01-3ba2-4d6b-a15d-e7048f78f22e",
        filename: "local.rds",
        filetype: "RDS",
        id: "8d4ccaf6-1ed7-4eed-b653-7206669209fa",
        s3_uri:
          "s3://corpora-data-prod/d4e69e01-3ba2-4d6b-a15d-e7048f78f22e/local.rds",
        type: "REMIX",
        updated_at: 1617323095.801812,
        user_submitted: true,
      },
    ],
    dataset_deployments: [
      {
        dataset_id: "d4e69e01-3ba2-4d6b-a15d-e7048f78f22e",
        url: "https://cellxgene.cziscience.com/e/d4e69e01-3ba2-4d6b-a15d-e7048f78f22e.cxg/",
      },
    ],
    development_stage: [
      {
        label: "40-year-old human stage",
        ontology_term_id: "HsapDv:0000134",
      },
      {
        label: "45-year-old human stage",
        ontology_term_id: "HsapDv:0000139",
      },
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "55-year-old human stage",
        ontology_term_id: "HsapDv:0000149",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "65-year-old human stage",
        ontology_term_id: "HsapDv:0000159",
      },
      {
        label: "70-year-old human stage",
        ontology_term_id: "HsapDv:0000164",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "d4e69e01-3ba2-4d6b-a15d-e7048f78f22e",
    is_valid: false,
    linked_genesets: [],
    name: "All \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    processing_status: {
      h5ad_status: "CONVERTED",
      cxg_status: "CONVERTED",
      rds_status: "CONVERTED",
      created_at: 1617319340.176877,
      dataset_id: "d4e69e01-3ba2-4d6b-a15d-e7048f78f22e",
      id: "02ad10b9-9580-4e04-a39d-d6ee7ccff69c",
      processing_status: "SUCCESS",
      updated_at: 1617323095.913868,
      upload_progress: 1.0,
      upload_status: "UPLOADED",
      validation_status: "VALID",
    },
    published: true,
    revision: 0,
    sex: ["female", "male"],
    tissue: [
      {
        label: "apex of heart",
        ontology_term_id: "UBERON:0002098",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "interventricular septum",
        ontology_term_id: "UBERON:0002094",
      },
      {
        label: "left cardiac atrium",
        ontology_term_id: "UBERON:0002079",
      },
      {
        label: "right cardiac atrium",
        ontology_term_id: "UBERON:0002078",
      },
    ],
    updated_at: 1617755780.2454,
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 40868,
    collection_id: "b52eb423-5d0d-4645-b217-e1c6d38b2e72",
    collection_visibility: "PUBLIC",
    created_at: 1617746898.168706,
    dataset_assets: [
      {
        created_at: 1617747713.055632,
        dataset_id: "ed852810-a003-4386-9846-1638362cee39",
        filename: "local.h5ad",
        filetype: "H5AD",
        id: "f46b664e-e415-4db9-8b89-f9bb18eb9ba5",
        s3_uri:
          "s3://corpora-data-prod/ed852810-a003-4386-9846-1638362cee39/local.h5ad",
        type: "REMIX",
        updated_at: 1617747713.055638,
        user_submitted: true,
      },
      {
        created_at: 1617748062.855596,
        dataset_id: "ed852810-a003-4386-9846-1638362cee39",
        filename: "local.rds",
        filetype: "RDS",
        id: "24960190-1eb9-4afa-ac26-d55e6981af4c",
        s3_uri:
          "s3://corpora-data-prod/ed852810-a003-4386-9846-1638362cee39/local.rds",
        type: "REMIX",
        updated_at: 1617748062.855601,
        user_submitted: true,
      },
    ],
    dataset_deployments: [
      {
        dataset_id: "ed852810-a003-4386-9846-1638362cee39",
        url: "https://cellxgene.cziscience.com/e/ed852810-a003-4386-9846-1638362cee39.cxg/",
      },
    ],
    development_stage: [
      {
        label: "40-year-old human stage",
        ontology_term_id: "HsapDv:0000134",
      },
      {
        label: "45-year-old human stage",
        ontology_term_id: "HsapDv:0000139",
      },
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "55-year-old human stage",
        ontology_term_id: "HsapDv:0000149",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "65-year-old human stage",
        ontology_term_id: "HsapDv:0000159",
      },
      {
        label: "70-year-old human stage",
        ontology_term_id: "HsapDv:0000164",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "ed852810-a003-4386-9846-1638362cee39",
    is_valid: false,
    linked_genesets: [],
    name: "Immune \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    processing_status: {
      h5ad_status: "CONVERTED",
      cxg_status: "CONVERTED",
      rds_status: "CONVERTED",
      created_at: 1617746898.171688,
      dataset_id: "ed852810-a003-4386-9846-1638362cee39",
      id: "18ab54b7-9c7d-4540-880b-e6e313cae131",
      processing_status: "SUCCESS",
      updated_at: 1617748062.928275,
      upload_progress: 1.0,
      upload_status: "UPLOADED",
      validation_status: "VALID",
    },
    published: true,
    revision: 0,
    sex: ["female", "male"],
    tissue: [
      {
        label: "apex of heart",
        ontology_term_id: "UBERON:0002098",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "interventricular septum",
        ontology_term_id: "UBERON:0002094",
      },
      {
        label: "left cardiac atrium",
        ontology_term_id: "UBERON:0002079",
      },
      {
        label: "right cardiac atrium",
        ontology_term_id: "UBERON:0002078",
      },
    ],
    updated_at: 1617755780.245403,
  },
  {
    assay: [
      {
        label: "10X 3' v2 sequencing",
        ontology_term_id: "EFO:0009899",
      },
      {
        label: "10x 3' v3 sequencing",
        ontology_term_id: "EFO:0009922",
      },
    ],
    cell_count: 3799,
    collection_id: "b52eb423-5d0d-4645-b217-e1c6d38b2e72",
    collection_visibility: "PUBLIC",
    created_at: 1617739319.906348,
    dataset_assets: [
      {
        created_at: 1617739451.659026,
        dataset_id: "f75f2ff4-2884-4c2d-b375-70de37a34507",
        filename: "local.h5ad",
        filetype: "H5AD",
        id: "323f4fca-f1a3-4921-a350-b2d33ca76df1",
        s3_uri:
          "s3://corpora-data-prod/f75f2ff4-2884-4c2d-b375-70de37a34507/local.h5ad",
        type: "REMIX",
        updated_at: 1617739451.659032,
        user_submitted: true,
      },
      {
        created_at: 1617739498.163832,
        dataset_id: "f75f2ff4-2884-4c2d-b375-70de37a34507",
        filename: "local.rds",
        filetype: "RDS",
        id: "a33fd442-e535-4e0b-a7e2-c26df804f9ab",
        s3_uri:
          "s3://corpora-data-prod/f75f2ff4-2884-4c2d-b375-70de37a34507/local.rds",
        type: "REMIX",
        updated_at: 1617739498.163839,
        user_submitted: true,
      },
    ],
    dataset_deployments: [
      {
        dataset_id: "f75f2ff4-2884-4c2d-b375-70de37a34507",
        url: "https://cellxgene.cziscience.com/e/f75f2ff4-2884-4c2d-b375-70de37a34507.cxg/",
      },
    ],
    development_stage: [
      {
        label: "40-year-old human stage",
        ontology_term_id: "HsapDv:0000134",
      },
      {
        label: "45-year-old human stage",
        ontology_term_id: "HsapDv:0000139",
      },
      {
        label: "50-year-old human stage",
        ontology_term_id: "HsapDv:0000144",
      },
      {
        label: "55-year-old human stage",
        ontology_term_id: "HsapDv:0000149",
      },
      {
        label: "60-year-old human stage",
        ontology_term_id: "HsapDv:0000154",
      },
      {
        label: "65-year-old human stage",
        ontology_term_id: "HsapDv:0000159",
      },
      {
        label: "70-year-old human stage",
        ontology_term_id: "HsapDv:0000164",
      },
    ],
    disease: [
      {
        label: "normal",
        ontology_term_id: "PATO:0000461",
      },
    ],
    ethnicity: [
      {
        label: "unknown",
        ontology_term_id: "",
      },
    ],
    id: "f75f2ff4-2884-4c2d-b375-70de37a34507",
    is_valid: false,
    linked_genesets: [],
    name: "Adipocytes \u2014 Cells of the adult human heart",
    organism: [
      {
        label: "Homo sapiens",
        ontology_term_id: "NCBITaxon:9606",
      },
    ],
    processing_status: {
      h5ad_status: "CONVERTED",
      cxg_status: "CONVERTED",
      rds_status: "CONVERTED",
      created_at: 1617739319.909173,
      dataset_id: "f75f2ff4-2884-4c2d-b375-70de37a34507",
      id: "6207defa-512b-4bbf-aba7-65911339aa6f",
      processing_status: "SUCCESS",
      updated_at: 1617739498.246613,
      upload_progress: 1.0,
      upload_status: "UPLOADED",
      validation_status: "VALID",
    },
    published: true,
    revision: 0,
    sex: ["female", "male"],
    tissue: [
      {
        label: "apex of heart",
        ontology_term_id: "UBERON:0002098",
      },
      {
        label: "heart left ventricle",
        ontology_term_id: "UBERON:0002084",
      },
      {
        label: "heart right ventricle",
        ontology_term_id: "UBERON:0002080",
      },
      {
        label: "interventricular septum",
        ontology_term_id: "UBERON:0002094",
      },
      {
        label: "left cardiac atrium",
        ontology_term_id: "UBERON:0002079",
      },
      {
        label: "right cardiac atrium",
        ontology_term_id: "UBERON:0002078",
      },
    ],
    updated_at: 1617755780.245404,
  },
];
export default datasets as unknown as Dataset[];
