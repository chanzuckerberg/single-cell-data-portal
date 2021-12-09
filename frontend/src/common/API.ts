export enum API {
  DATASET = "/dp/v1/datasets/{dataset_uuid}",
  DATASET_ASSETS = "/dp/v1/datasets/{dataset_uuid}/assets",
  DATASET_ASSET_DOWNLOAD_LINK = "/dp/v1/datasets/{dataset_uuid}/asset/{asset_uuid}",
  DATASET_STATUS = "/dp/v1/datasets/{dataset_uuid}/status",
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
}
