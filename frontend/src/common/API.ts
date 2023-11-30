export enum API {
  CURATOR_AUTH_KEY = "/dp/v1/auth/key",
  DATASET = "/dp/v1/datasets/{dataset_id}",
  DATASET_ASSETS = "/dp/v1/datasets/{dataset_id}/assets",
  DATASET_ASSET_DOWNLOAD_LINK = "/dp/v1/datasets/{dataset_id}/asset/{asset_id}",
  DATASET_STATUS = "/dp/v1/datasets/{dataset_id}/status",
  DATASETS_INDEX = "/dp/v1/datasets/index", // Filter-specific endpoint
  COLLECTIONS_INDEX = "/dp/v1/collections/index", // Filter-specific endpoint
  COLLECTION = "/dp/v1/collections/{id}",
  COLLECTION_UPLOAD_LINKS = "/dp/v1/collections/{id}/upload-links",
  COLLECTION_PUBLISH = "/dp/v1/collections/{id}/publish",
  CREATE_COLLECTION = "/dp/v1/collections",
  LOG_IN = "/dp/v1/login",
  LOG_OUT = "/dp/v1/logout",
  USER_COLLECTIONS_INDEX = "/dp/v1/user-collections/index", // Filter-specific endpoint
  USER_DATASETS_INDEX = "/dp/v1/user-datasets/index", // Filter-specific endpoint
  USER_INFO = "/dp/v1/userinfo",
  WMG_PRIMARY_FILTER_DIMENSIONS = "/wmg/v1/primary_filter_dimensions",
  WMG_QUERY = "/wmg/v1/query",
  WMG_FILTERS_QUERY = "/wmg/v1/filters",
  WMG_MARKER_GENES = "/wmg/v2/markers",
  WMG_GENE_INFO = "/gene_info/v1/gene_info",
  CENSUS_DIRECTORY = "/census/v1/directory", // TODO: replace with actual endpoint
}
