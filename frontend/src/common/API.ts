export enum API {
  DATASET_ASSET_DOWNLOAD_LINK = "/dp/v1/dataset/{dataset_uuid}/asset/{asset_uuid}",
  COLLECTIONS = "/dp/v1/collections",
  COLLECTION = "/dp/v1/collections/{id}",
  COLLECTION_UPLOAD_LINKS = "/dp/v1/collections/{id}/upload-links",
  CREATE_COLLECTION = "/dp/v1/collections",
  LOG_IN = "/dp/v1/login",
  LOG_OUT = "/dp/v1/logout",
  USER_INFO = "/dp/v1/userinfo",
}
