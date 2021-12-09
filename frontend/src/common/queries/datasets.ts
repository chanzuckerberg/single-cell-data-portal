import { QueryResult, useMutation, useQuery, useQueryCache } from "react-query";
import { API_URL } from "src/configs/configs";
import { API } from "../API";
import {
  DatasetMetadata,
  DatasetUploadStatus,
  VISIBILITY_TYPE,
} from "../entities";
import { apiTemplateToUrl } from "../utils/apiTemplateToUrl";
import { USE_COLLECTION } from "./collections";
import { DEFAULT_FETCH_OPTIONS, DELETE_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";

export const USE_DATASET_STATUS = {
  entities: [ENTITIES.DATASET_STATUS],
  id: "dataset",
};

async function fetchDatasetStatus(
  _: unknown,
  dataset_uuid: string
): Promise<DatasetUploadStatus> {
  const url = apiTemplateToUrl(API_URL + API.DATASET_STATUS, { dataset_uuid });

  return await (await fetch(url, DEFAULT_FETCH_OPTIONS)).json();
}

const REFETCH_INTERVAL_MS = 10 * 1000;

export function useDatasetStatus(dataset_uuid: string, shouldFetch: boolean) {
  return useQuery<DatasetUploadStatus>(
    [USE_DATASET_STATUS, dataset_uuid],
    fetchDatasetStatus,
    { enabled: shouldFetch, refetchInterval: REFETCH_INTERVAL_MS }
  );
}

export const USE_DELETE_DATASET = {
  entities: [ENTITIES.DATASET],
  id: "dataset",
};

async function deleteDataset(dataset_uuid = ""): Promise<DatasetUploadStatus> {
  if (!dataset_uuid) throw new Error("No dataset id provided");

  const url = apiTemplateToUrl(API_URL + API.DATASET, { dataset_uuid });
  const response = await fetch(url, DELETE_FETCH_OPTIONS);

  if (response.ok) return await response.json();

  throw Error(response.statusText);
}

export function useDeleteDataset(collection_uuid = "") {
  if (!collection_uuid) {
    throw new Error("No collection id given");
  }

  const queryCache = useQueryCache();

  return useMutation(deleteDataset, {
    onSuccess: (uploadStatus: DatasetUploadStatus) => {
      queryCache.invalidateQueries([
        USE_COLLECTION,
        collection_uuid,
        VISIBILITY_TYPE.PRIVATE,
      ]);

      queryCache.cancelQueries([USE_DATASET_STATUS, uploadStatus.dataset_id]);

      queryCache.setQueryData(
        [USE_DATASET_STATUS, uploadStatus.dataset_id],
        uploadStatus
      );
    },
  });
}

/**
 * Query key for /dataset-metadata/id
 */
export const USE_DATASET_METADATA = {
  entities: [ENTITIES.DATASET],
  id: "datasetMetadata",
};

/**
 * Cache-enabled hook for fetching metadata for the dataset with the given explorer URL.
 * @param explorerUrl - Explorer URL of a dataset.
 * @param enabled - True if fetch can be invoked.
 * @returns Dataset metadata.
 */
export function useFetchDatasetMetadata(
  explorerUrl: string,
  enabled: boolean
): QueryResult<DatasetMetadata> {
  return useQuery<DatasetMetadata>(
    [USE_DATASET_METADATA, explorerUrl],
    () => fetchDatasetMetadata(explorerUrl),
    { enabled }
  );
}

/**
 * Fetch metadata for the dataset with the given explorer URL.
 * @param explorerUrl - Explorer URL of a dataset.
 * @returns Promise that resolves to dataset metadata.
 */
async function fetchDatasetMetadata(
  explorerUrl: string
): Promise<DatasetMetadata> {
  const explorerPath = new URL(explorerUrl).pathname;
  const { metadata } = await (
    await fetch(
      apiTemplateToUrl(API_URL + API.DATASET_METADATA, { explorerPath })
    )
  ).json();
  return metadata;
}
