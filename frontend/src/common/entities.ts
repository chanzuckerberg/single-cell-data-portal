export enum COLLECTION_LINK_TYPE {
  DOI = "DOI",
  RAW_DATA = "RAW_DATA",
  PROTOCOL = "PROTOCOL",
  LAB_WEBSITE = "LAB_WEBSITE",
  OTHER = "OTHER",
}

export const COLLECTION_LINK_TYPE_OPTIONS = {
  [COLLECTION_LINK_TYPE.DOI]: { text: "DOI", value: COLLECTION_LINK_TYPE.DOI },
  [COLLECTION_LINK_TYPE.RAW_DATA]: {
    text: "Raw Data",
    value: COLLECTION_LINK_TYPE.RAW_DATA,
  },
  [COLLECTION_LINK_TYPE.PROTOCOL]: {
    text: "Protocol",
    value: COLLECTION_LINK_TYPE.PROTOCOL,
  },
  [COLLECTION_LINK_TYPE.LAB_WEBSITE]: {
    text: "Lab Website",
    value: COLLECTION_LINK_TYPE.LAB_WEBSITE,
  },
  [COLLECTION_LINK_TYPE.OTHER]: {
    text: "Other",
    value: COLLECTION_LINK_TYPE.OTHER,
  },
};

export interface Link {
  link_name: string;
  link_url: string;
  link_type: COLLECTION_LINK_TYPE;
}

export enum ACCESS_TYPE {
  READ = "READ",
  WRITE = "WRITE",
}

export enum VISIBILITY_TYPE {
  PUBLIC = "PUBLIC",
  PRIVATE = "PRIVATE",
}

export interface Collection {
  access_type: ACCESS_TYPE;
  contact_email: string;
  contact_name: string;
  description: string;
  id: string;
  organs: string[];
  name: string;
  owner: string;
  visibility: VISIBILITY_TYPE;
  datasets: Dataset[];
  links: Link[];
  data_submission_policy_version: string;
  obfuscated_uuid: string;
  created_at: number;
  updated_at: number;
}

export interface Dataset {
  id: string;
  // assay: string;
  // assay_ontology: string;
  // tissue: string;
  // tissue_ontology: string;
  // disease_state: string;
  // disease_state_ontology: string;
  // sex: string;
  // ethnicity: string;
  // ethnicity_ontology: string;
  // organism: string;
  // organism_ontology: string;
  name: string;
  // source_data_location: string;
  // revision: number;
  dataset_deployments: DatasetDeployment[];
  dataset_assets: DatasetAsset[];
  // contributors: Contributor[];
  // preprint_doi: DOI;
  // publication_doi: DOI;
}

export enum DATASET_ASSET_FORMAT {
  H5AD = "H5AD",
  RDS = "RDS",
  LOOM = "LOOM",
  CXG = "CXG",
}

export enum DATASET_ASSET_TYPE {
  REMIX = "REMIX",
  ORIGINAL = "ORIGINAL",
}

export interface DatasetAsset {
  id: string;
  dataset_id: string;
  filetype: DATASET_ASSET_FORMAT;
  filename: string;
  type: DATASET_ASSET_TYPE;
  s3_uri: string;
}

export interface DatasetDeployment {
  id: string;
  environment: string;
  url: string;
}

// interface Contributor {
//   id: string;
//   name: string;
//   email: string;
//   institution: string;
// }

// interface DOI {
//   id: string;
//   title: string;
//   data: string;
//   journal: string;
// }

export interface User {
  picture: string;
}
