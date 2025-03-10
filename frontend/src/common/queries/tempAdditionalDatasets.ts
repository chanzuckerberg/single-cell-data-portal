/* eslint-disable sonarjs/no-duplicate-string */

import { TISSUE_TYPE } from "src/components/common/Filter/common/entities";
import { IS_PRIMARY_DATA } from "src/common/entities";

export const newWormDataset = {
  assay: [
    {
      label: "10x 3' v2",
      ontology_term_id: "EFO:0009899",
    },
  ],
  cell_count: 8362,
  cell_type: [
    {
      label: "HOso Cell",
      ontology_term_id: "WBbt:0004759",
    },
  ],
  cell_type_ancestors: ["WBbt:0006930", "WBbt:0005750", "WBbt:0005762"],
  collection_id: "89f3e0a4-bbe5-4861-9add-71526aa62afd",
  created_at: 1728319480.965303,
  dataset_assets: [
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "raw.h5ad",
      filetype: "RAW_H5AD",
      id: "06ab3c80-2d8a-4f70-b605-63725c9bb52f",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/raw.h5ad",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "local.h5ad",
      filetype: "H5AD",
      id: "037f6443-6e21-4928-82f9-e0a9a79a04c0",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/local.h5ad",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "local.rds",
      filetype: "RDS",
      id: "56a5b2ba-a612-4c38-8db0-1d79b58f627a",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/local.rds",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "",
      filetype: "CXG",
      id: "ce439a58-e16b-4da2-b23b-5f073639f509",
      s3_uri:
        "s3://hosted-cellxgene-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596.cxg/",
      updated_at: 0,
      user_submitted: true,
    },
  ],
  dataset_deployments: [
    {
      url: "https://cellxgene.dev.single-cell.czi.technology/e/cd77258f-b08b-4c89-b93f-6e6f146b1a4d.cxg/",
    },
  ],
  development_stage: [
    {
      label: "10-days post-L4 adult hermaphrodite Ce",
      ontology_term_id: "WBls:0000674",
    },
    {
      label: "11-days post-L4 adult hermaphrodite Ce",
      ontology_term_id: "WBls:0000797",
    },
    {
      label: "12-days post-L4 adult hermaphrodite Ce",
      ontology_term_id: "WBls:0000798",
    },
  ],
  development_stage_ancestors: ["WBls:0000804", "WBls:0000816"],
  disease: [
    {
      label: "Alzheimer disease",
      ontology_term_id: "MONDO:0004975",
    },
    {
      label: "normal",
      ontology_term_id: "PATO:0000461",
    },
  ],
  donor_id: ["3", "1", "2", "5", "6", "7", "4", "8", "9", "10"],
  explorer_url:
    "https://cellxgene.dev.single-cell.czi.technology/e/cd77258f-b08b-4c89-b93f-6e6f146b1a4d.cxg/",
  id: "000198ac-27c7-4b9f-9fd4-6362a5b91594",
  is_primary_data: IS_PRIMARY_DATA.SECONDARY,
  is_valid: true,
  mean_genes_per_cell: 2165.4574264530015,
  name: "FAKE FAKE FAKE - Worm HOso Cell",
  organism: [
    {
      label: "Caenorhabditis elegans",
      ontology_term_id: "NCBITaxon:6239",
    },
  ],
  primary_cell_count: 0,
  processing_status: {
    created_at: 0,
    cxg_status: "UPLOADED",
    dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
    h5ad_status: "UPLOADED",
    id: "NA",
    processing_status: "SUCCESS",
    rds_status: "UPLOADED",
    updated_at: 0,
    upload_progress: 1,
    upload_status: "UPLOADED",
    validation_status: "VALID",
  },
  published: true,
  published_at: 1605879580.522766,
  revision: 0,
  schema_version: "5.2.0",
  self_reported_ethnicity: [
    {
      label: "unknown",
      ontology_term_id: "unknown",
    },
  ],
  sex: [
    {
      label: "male",
      ontology_term_id: "PATO:0000384",
    },
  ],
  suspension_type: ["nucleus"],
  tissue: [
    {
      label: "hook sensillum",
      ontology_term_id: "WBbt:0006930",
      tissue_type: TISSUE_TYPE.TISSUE,
    },
  ],
  tissue_ancestors: ["WBbt:0006929", "WBbt:0003760"],
  tombstone: false,
  updated_at: 1728319480.965303,
};

export const newFruitFlyDataset = {
  assay: [
    {
      label: "10x 3' v2",
      ontology_term_id: "EFO:0009899",
    },
  ],
  cell_count: 8362,
  cell_type: [
    {
      label: "germline stem cell",
      ontology_term_id: "FBbt:00004861",
    },
  ],
  cell_type_ancestors: ["FBbt:00004860"],
  collection_id: "89f3e0a4-bbe5-4861-9add-71526aa62afd",
  created_at: 1728319480.965303,
  dataset_assets: [
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "raw.h5ad",
      filetype: "RAW_H5AD",
      id: "06ab3c80-2d8a-4f70-b605-63725c9bb52f",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/raw.h5ad",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "local.h5ad",
      filetype: "H5AD",
      id: "037f6443-6e21-4928-82f9-e0a9a79a04c0",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/local.h5ad",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "local.rds",
      filetype: "RDS",
      id: "56a5b2ba-a612-4c38-8db0-1d79b58f627a",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/local.rds",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "",
      filetype: "CXG",
      id: "ce439a58-e16b-4da2-b23b-5f073639f509",
      s3_uri:
        "s3://hosted-cellxgene-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596.cxg/",
      updated_at: 0,
      user_submitted: true,
    },
  ],
  dataset_deployments: [
    {
      url: "https://cellxgene.dev.single-cell.czi.technology/e/cd77258f-b08b-4c89-b93f-6e6f146b1a4d.cxg/",
    },
  ],
  development_stage: [
    {
      label: "day 11 of adulthood",
      ontology_term_id: "FBdv:00007086",
    },
  ],
  development_stage_ancestors: ["FBdv:00007002"],
  disease: [
    {
      label: "Alzheimer disease",
      ontology_term_id: "MONDO:0004975",
    },
    {
      label: "normal",
      ontology_term_id: "PATO:0000461",
    },
  ],
  donor_id: ["3", "1", "2", "5", "6", "7", "4", "8", "9", "10"],
  explorer_url:
    "https://cellxgene.dev.single-cell.czi.technology/e/cd77258f-b08b-4c89-b93f-6e6f146b1a4d.cxg/",
  id: "000198ac-27c7-4b9f-9fd4-6362a5b91597",
  is_primary_data: IS_PRIMARY_DATA.SECONDARY,
  is_valid: true,
  mean_genes_per_cell: 2165.4574264530015,
  name: "FAKE FAKE FAKE - new Fruit Fly - Drosophila melanogaster",
  organism: [
    {
      label: "Drosophila melanogaster",
      ontology_term_id: "NCBITaxon:7227",
    },
  ],
  primary_cell_count: 0,
  processing_status: {
    created_at: 0,
    cxg_status: "UPLOADED",
    dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
    h5ad_status: "UPLOADED",
    id: "NA",
    processing_status: "SUCCESS",
    rds_status: "UPLOADED",
    updated_at: 0,
    upload_progress: 1,
    upload_status: "UPLOADED",
    validation_status: "VALID",
  },
  published: true,
  published_at: 1605879580.522766,
  revision: 0,
  schema_version: "5.2.0",
  self_reported_ethnicity: [
    {
      label: "unknown",
      ontology_term_id: "unknown",
    },
  ],
  sex: [
    {
      label: "female",
      ontology_term_id: "PATO:0000383",
    },
  ],
  suspension_type: ["nucleus"],
  tissue: [
    {
      label: "NB6-2 lineage neuron",
      ontology_term_id: "FBbt:00048943",
      tissue_type: TISSUE_TYPE.TISSUE,
    },
  ],
  tissue_ancestors: ["FBbt:00025990", "FBbt:00000000"],
  tombstone: false,
  updated_at: 1728319480.965303,
};

export const newZebraFishDataset = {
  assay: [
    {
      label: "10x 3' v2",
      ontology_term_id: "EFO:0009899",
    },
  ],
  cell_count: 8362,
  cell_type: [
    {
      label: "dental papilla cell",
      ontology_term_id: "ZFA:0009173",
    },
  ],
  cell_type_ancestors: ["ZFA:0005140"],
  collection_id: "89f3e0a4-bbe5-4861-9add-71526aa62afd",
  created_at: 1728319480.965303,
  dataset_assets: [
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "raw.h5ad",
      filetype: "RAW_H5AD",
      id: "06ab3c80-2d8a-4f70-b605-63725c9bb52f",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/raw.h5ad",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "local.h5ad",
      filetype: "H5AD",
      id: "037f6443-6e21-4928-82f9-e0a9a79a04c0",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/local.h5ad",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "local.rds",
      filetype: "RDS",
      id: "56a5b2ba-a612-4c38-8db0-1d79b58f627a",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/local.rds",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "",
      filetype: "CXG",
      id: "ce439a58-e16b-4da2-b23b-5f073639f509",
      s3_uri:
        "s3://hosted-cellxgene-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596.cxg/",
      updated_at: 0,
      user_submitted: true,
    },
  ],
  dataset_deployments: [
    {
      url: "https://cellxgene.dev.single-cell.czi.technology/e/cd77258f-b08b-4c89-b93f-6e6f146b1a4d.cxg/",
    },
  ],
  development_stage: [
    {
      label: "Hatching:Long-pec",
      ontology_term_id: "ZFS:0000033",
    },
  ],
  development_stage_ancestors: ["ZFS:0000007"],
  disease: [
    {
      label: "Alzheimer disease",
      ontology_term_id: "MONDO:0004975",
    },
    {
      label: "normal",
      ontology_term_id: "PATO:0000461",
    },
  ],
  donor_id: ["3", "1", "2", "5", "6", "7", "4", "8", "9", "10"],
  explorer_url:
    "https://cellxgene.dev.single-cell.czi.technology/e/cd77258f-b08b-4c89-b93f-6e6f146b1a4d.cxg/",
  id: "000198ac-27c7-4b9f-9fd4-6362a5b91591",
  is_primary_data: IS_PRIMARY_DATA.SECONDARY,
  is_valid: true,
  mean_genes_per_cell: 2165.4574264530015,
  name: "FAKE FAKE FAKE - Zebrafish - Danio rerio",
  organism: [
    {
      label: "Danio rerio",
      ontology_term_id: "NCBITaxon:7955",
    },
  ],
  primary_cell_count: 0,
  processing_status: {
    created_at: 0,
    cxg_status: "UPLOADED",
    dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
    h5ad_status: "UPLOADED",
    id: "NA",
    processing_status: "SUCCESS",
    rds_status: "UPLOADED",
    updated_at: 0,
    upload_progress: 1,
    upload_status: "UPLOADED",
    validation_status: "VALID",
  },
  published: true,
  published_at: 1605879580.522766,
  revision: 0,
  schema_version: "5.2.0",
  self_reported_ethnicity: [
    {
      label: "unknown",
      ontology_term_id: "unknown",
    },
  ],
  sex: [
    {
      label: "male",
      ontology_term_id: "PATO:0000384",
    },
  ],
  suspension_type: ["nucleus"],
  tissue: [
    {
      label: "dental papilla",
      ontology_term_id: "ZFA:0005140",
      tissue_type: TISSUE_TYPE.TISSUE,
    },
  ],
  tissue_ancestors: ["ZFA:0001477", "ZFA:0000037"],
  tombstone: false,
  updated_at: 1728319480.965303,
};

export const newFruitFlyDataset2 = {
  assay: [
    {
      label: "10x 3' v2",
      ontology_term_id: "EFO:0009899",
    },
  ],
  cell_count: 8362,
  cell_type: [
    {
      label: "larval abdominal fat cell",
      ontology_term_id: "FBbt:00049952",
    },
  ],
  cell_type_ancestors: ["FBbt:00004860"],
  collection_id: "89f3e0a4-bbe5-4861-9add-71526aa62afd",
  created_at: 1728319480.965303,
  dataset_assets: [
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "raw.h5ad",
      filetype: "RAW_H5AD",
      id: "06ab3c80-2d8a-4f70-b605-63725c9bb52f",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/raw.h5ad",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "local.h5ad",
      filetype: "H5AD",
      id: "037f6443-6e21-4928-82f9-e0a9a79a04c0",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/local.h5ad",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "local.rds",
      filetype: "RDS",
      id: "56a5b2ba-a612-4c38-8db0-1d79b58f627a",
      s3_uri:
        "s3://corpora-data-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596/local.rds",
      updated_at: 0,
      user_submitted: true,
    },
    {
      created_at: 0,
      dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
      filename: "",
      filetype: "CXG",
      id: "ce439a58-e16b-4da2-b23b-5f073639f509",
      s3_uri:
        "s3://hosted-cellxgene-dev/000198ac-27c7-4b9f-9fd4-6362a5b91596.cxg/",
      updated_at: 0,
      user_submitted: true,
    },
  ],
  dataset_deployments: [
    {
      url: "https://cellxgene.dev.single-cell.czi.technology/e/cd77258f-b08b-4c89-b93f-6e6f146b1a4d.cxg/",
    },
  ],
  development_stage: [
    {
      label: "Adult",
      ontology_term_id: "FBdv:00007002",
    },
  ],
  development_stage_ancestors: [
    "FBdv:00005342",
    "FBdv:00005372",
    "FBdv:00005336",
  ],
  disease: [
    {
      label: "normal",
      ontology_term_id: "PATO:0000461",
    },
  ],
  donor_id: ["3", "1", "2", "5", "6", "7", "4", "8", "9", "10"],
  explorer_url:
    "https://cellxgene.dev.single-cell.czi.technology/e/cd77258f-b08b-4c89-b93f-6e6f146b1a4d.cxg/",
  id: "000198ac-27c7-4b9f-9fd4-6362a5b91597",
  is_primary_data: IS_PRIMARY_DATA.SECONDARY,
  is_valid: true,
  mean_genes_per_cell: 2165.4574264530015,
  name: "FAKE - new Fruit Fly 2  - Drosophila melanogaster",
  organism: [
    {
      label: "Drosophila melanogaster",
      ontology_term_id: "NCBITaxon:7227",
    },
  ],
  primary_cell_count: 0,
  processing_status: {
    created_at: 0,
    cxg_status: "UPLOADED",
    dataset_id: "000198ac-27c7-4b9f-9fd4-6362a5b91596",
    h5ad_status: "UPLOADED",
    id: "NA",
    processing_status: "SUCCESS",
    rds_status: "UPLOADED",
    updated_at: 0,
    upload_progress: 1,
    upload_status: "UPLOADED",
    validation_status: "VALID",
  },
  published: true,
  published_at: 1605879580.522766,
  revision: 0,
  schema_version: "5.2.0",
  self_reported_ethnicity: [
    {
      label: "unknown",
      ontology_term_id: "unknown",
    },
  ],
  sex: [
    {
      label: "male",
      ontology_term_id: "PATO:0000384",
    },
  ],
  suspension_type: ["nucleus"],
  tissue: [
    {
      label: "somatic cell of fat body/gonad primordium",
      ontology_term_id: "FBbt:00059345",
      tissue_type: TISSUE_TYPE.TISSUE,
    },
  ],
  tissue_ancestors: ["FBbt:00005520", "FBbt:00006001", "FBbt:00026000"],
  tombstone: false,
  updated_at: 1728319480.965303,
};
