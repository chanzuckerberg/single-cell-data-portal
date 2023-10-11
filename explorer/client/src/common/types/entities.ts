// If a globally shared type or interface doesn't have a clear owner, put it here

/**
 * Complete list of standard cellxgene column names.
 * https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md#obs-cell-metadata
 */
export const STANDARD_CATEGORY_NAMES = [
  "assay",
  "assay_ontology_term_id",
  "cell_type",
  "cell_type_ontology_term_id",
  "development_stage",
  "development_stage_ontology_term_id",
  "disease",
  "disease_ontology_term_id",
  "self_reported_ethnicity",
  "self_reported_ethnicity_ontology_term_id",
  "is_primary_data",
  "organism",
  "organism_ontology_term_id",
  "sex_ontology_term_id",
  "sex",
  "tissue_ontology_term_id",
  "tissue",
  "suspension_type",
  "donor_id",
];

/**
 * Author of publication associated with a collection, populated from Crossref as part of collection publication
 * metadata.
 */
export interface Author {
  family: string;
  given: string;
}

/**
 * Consortium of publication associated with a collection, populated from Crossref as part of collection publication
 * metadata.
 */
export interface Consortium {
  name: string;
}

/**
 * Collection and dataset information for the selected dataset.
 */
export interface DatasetMetadata {
  collection_contact_email: string;
  collection_contact_name: string;
  collection_datasets: Dataset[];
  collection_description: string;
  collection_links: Link[];
  collection_name: string;
  collection_publisher_metadata: PublisherMetadata;
  collection_url: string;
  dataset_id: string;
  dataset_name: string;
  s3_URI: S3URI;
}

export type S3URI = string;

/**
 * Corpora encoding and schema versioning.
 */
interface CorporaVersion {
  corpora_encoding_version: string;
  corpora_schema_version: string;
}

/**
 * Basic dataset model.
 */
export interface Dataset {
  id: string;
  cell_count: number;
  dataset_deployments: DatasetDeployment[];
  name: string;
}

/**
 * Deployment location of dataset.
 */
export interface DatasetDeployment {
  url: string;
}

/**
 * Flags informing garbage collection-related logic.
 */
export interface GCHints {
  isHot: boolean;
}

/**
 * Rough Data Portal values returned in config.
 */
export interface DataPortalProps {
  version: CorporaVersion;
  title: string;
  contributors: { name: string; institution: string }[];
  layer_descriptions: { X: string };
  organism: string;
  organism_ontology_term_id: string;
  project_name: string;
  project_description: string;
  project_links: { link_name: string; link_type: string; link_url: string }[];
  default_embedding: string;
  preprint_doi: string;
  publication_doi: string;
  schema_version: string;
}

/**
 * Type and location of a related resource.
 */
export interface Link {
  link_name: string;
  link_type: string;
  link_url: string;
}

/**
 * Collection publication metadata, populated from Crossref by collection publication DOI.
 */
export interface PublisherMetadata {
  authors: (Author | Consortium)[];
  journal: string;
  published_day: number;
  published_month: number;
  published_year: number;
}
