import { useMutation, useQuery, useQueryCache } from "react-query";
import { API_URL } from "src/configs/configs";
import { API } from "../API";
import { DatasetUploadStatus, VISIBILITY_TYPE } from "../entities";
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

async function deleteDataset(dataset_uuid: string): Promise<boolean> {
  const url = apiTemplateToUrl(API_URL + API.DATASET, { dataset_uuid });
  return await (await fetch(url, DELETE_FETCH_OPTIONS)).json();
}

export function useDeleteDataset(
  dataset_uuid: string,
  collection_uuid: string
) {
  const queryCache = useQueryCache();

  return useMutation(deleteDataset, {
    onSuccess: () => {
      queryCache.invalidateQueries([
        USE_COLLECTION,
        collection_uuid,
        VISIBILITY_TYPE.PRIVATE,
      ]);
      queryCache.cancelQueries([USE_DATASET_STATUS, dataset_uuid]);
    },
  });
}
