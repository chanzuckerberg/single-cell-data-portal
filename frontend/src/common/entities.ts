import { CONSORTIA } from "src/components/CreateCollectionModal/components/Content/common/constants";

export enum COLLECTION_LINK_TYPE {
  DOI = "DOI",
  RAW_DATA = "RAW_DATA",
  PROTOCOL = "PROTOCOL",
  LAB_WEBSITE = "LAB_WEBSITE",
  OTHER = "OTHER",
  DATA_SOURCE = "DATA_SOURCE",
}

/**
 * Options backing link type menu, used when associating links with collection.
 */
export const COLLECTION_LINK_TYPE_OPTIONS = {
  [COLLECTION_LINK_TYPE.DOI]: {
    text: "Publication DOI",
    value: COLLECTION_LINK_TYPE.DOI,
  },
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
  [COLLECTION_LINK_TYPE.DATA_SOURCE]: {
    text: "Data Source",
    value: COLLECTION_LINK_TYPE.DATA_SOURCE,
  },
};

/**
 * Collection status, used to determine whether a collection is private or public, or in revision.
 */
export enum COLLECTION_STATUS {
  PRIVATE = "Private",
  PUBLISHED = "Published",
  REVISION = "Revision",
}

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

export enum X_APPROXIMATE_DISTRIBUTION {
  COUNT = "count",
  NORMAL = "normal",
}

export enum IS_PRIMARY_DATA {
  PRIMARY = "PRIMARY",
  SECONDARY = "SECONDARY",
  BOTH = "BOTH",
}

/**
 * Author of publication associated with a collection, populated from Crossref as part of collection publication
 * metadata.
 */
export interface Author {
  family: string;
  given: string;
}

export interface Collection {
  access_type: ACCESS_TYPE;
  consortia: CONSORTIA[];
  contact_email: string;
  contact_name: string;
  description: string;
  id: string;
  organs: string[];
  name: string;
  owner: string;
  visibility: VISIBILITY_TYPE;
  datasets: Map<Dataset["id"], Dataset>;
  links: Link[];
  data_submission_policy_version: string;
  created_at: number;
  updated_at: number;
  publisher_metadata: PublisherMetadata;
  revision_diff: boolean;
  summaryCitation?: string;
  tombstone?: boolean;
  revising_in?: Collection["id"];
  revision_of?: Collection["id"];
}

/**
 * Consortium of publication associated with a collection, populated from Crossref as part of collection publication
 * metadata.
 */
export interface Consortium {
  name: string;
}

export type Ontology = {
  label: string;
  ontology_term_id: string;
};

export interface Dataset {
  id: string;
  assay: Ontology[];
  tissue: Ontology[];
  disease: Ontology[];
  cell_count: number | null;
  // sex: string;
  self_reported_ethnicity: Ontology;
  organism: Ontology[];
  name: string;
  cell_type: Ontology[];
  is_primary_data: IS_PRIMARY_DATA;
  x_approximate_distribution?: X_APPROXIMATE_DISTRIBUTION;
  schema_version: string;
  // source_data_location: string;
  // revision: number;
  dataset_deployments: DatasetDeployment[];
  dataset_assets: DatasetAsset[];
  processing_status: DatasetUploadStatus;
  collection_id: Collection["id"];
  // contributors: Contributor[];
  // preprint_doi: DOI;
  // publication_doi: DOI;
  created_at: number;
  original_id?: string;
  published?: boolean;
  published_at?: number;
  tombstone?: boolean;
  updated?: boolean;
  collection_visibility: Collection["visibility"];
}

export enum DATASET_ASSET_FORMAT {
  H5AD = "H5AD",
  RDS = "RDS",
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

export enum UPLOAD_STATUS {
  WAITING = "WAITING",
  UPLOADING = "UPLOADING",
  UPLOADED = "UPLOADED",
  FAILED = "FAILED",
  CANCEL_PENDING = "CANCEL_PENDING",
  PENDING = "PENDING",
  CANCELED = "CANCELED",
  NA = "NA",
}

export enum VALIDATION_STATUS {
  VALIDATING = "VALIDATING",
  VALID = "VALID",
  INVALID = "INVALID",
  NA = "NA",
}

export enum CONVERSION_STATUS {
  CONVERTING = "CONVERTING",
  CONVERTED = "CONVERTED",
  FAILED = "FAILED",
  SKIPPED = "SKIPPED",
  NA = "NA",
}

export enum PROCESSING_STATUS {
  PENDING = "PENDING",
  SUCCESS = "SUCCESS",
  FAILURE = "FAILURE",
}

export interface DatasetUploadStatus {
  dataset_id: string;
  upload_status: UPLOAD_STATUS;
  processing_status: PROCESSING_STATUS;
  upload_message: string;
  upload_progress: number;
  validation_status: VALIDATION_STATUS;
  validation_message: string;
  h5ad_status: CONVERSION_STATUS;
  cxg_status: CONVERSION_STATUS;
  rds_status: CONVERSION_STATUS;
}

export interface GeneSet {
  name: string;
  id: string;
  description: string;
  linked_datasets: Dataset["id"][];
}

/**
 * Collection publication metadata, populated from Crossref by collection publication DOI.
 */
export interface PublisherMetadata {
  authors: (Author | Consortium)[];
  journal: string;
  published_at: number;
  published_day: number;
  published_month: number;
  published_year: number;
}
