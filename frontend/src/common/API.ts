export enum API {
  CURATOR_AUTH_KEY = "/dp/v1/auth/key",
  DATASET = "/dp/v1/datasets/{dataset_id}",
  DATASET_ASSETS = "/dp/v1/datasets/{dataset_id}/assets",
  DATASET_ASSET_DOWNLOAD_LINK = "/dp/v1/datasets/{dataset_id}/asset/{asset_id}",
  DATASET_STATUS = "/dp/v1/datasets/{dataset_id}/status",
  DATASETS_INDEX = "/dp/v1/datasets/index", // Filter-specific endpoint
  COLLECTIONS = "/dp/v1/collections",
  COLLECTIONS_INDEX = "/dp/v1/collections/index", // Filter-specific endpoint
  COLLECTION = "/dp/v1/collections/{id}",
  COLLECTION_UPLOAD_LINKS = "/dp/v1/collections/{id}/upload-links",
  COLLECTION_PUBLISH = "/dp/v1/collections/{id}/publish",
  CREATE_COLLECTION = "/dp/v1/collections",
  LOG_IN = "/dp/v1/login",
  LOG_OUT = "/dp/v1/logout",
  USER_INFO = "/dp/v1/userinfo",
  WMG_PRIMARY_FILTER_DIMENSIONS = "/wmg/v1/primary_filter_dimensions",
  WMG_QUERY = "/wmg/v1/query",
}
