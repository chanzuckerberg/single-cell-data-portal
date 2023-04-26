import {
  useMutation,
  UseMutationResult,
  useQuery,
  useQueryClient,
  UseQueryResult,
} from "react-query";
import { API_URL } from "src/configs/configs";
import { API } from "../API";
import {
  Collection,
  Dataset,
  DatasetAsset,
  DatasetUploadStatus,
  PROCESSING_STATUS,
} from "../entities";
import { apiTemplateToUrl } from "../utils/apiTemplateToUrl";
import { USE_COLLECTION } from "./collections";
import { DEFAULT_FETCH_OPTIONS, DELETE_FETCH_OPTIONS } from "./common";
import { ENTITIES } from "./entities";
import { USE_DATASETS_INDEX } from "src/common/queries/filter";

/**
 * Cached query matching the refetch predicate, that are not being rendered, will be invalidated and refetched
 * in the background.
 */
const DEFAULT_BACKGROUND_REFETCH = {
  refetchInactive: true,
};

export const USE_DATASET_STATUS = {
  entities: [ENTITIES.DATASET_STATUS],
  id: "dataset",
};

async function fetchDatasetStatus(
  dataset_id: string
): Promise<DatasetUploadStatus> {
  const url = apiTemplateToUrl(API_URL + API.DATASET_STATUS, { dataset_id });

  return await (await fetch(url, DEFAULT_FETCH_OPTIONS)).json();
}

const REFETCH_INTERVAL_MS = 10 * 1000;

export function useDatasetStatus(
  dataset_id: string,
  collectionId: Collection["id"],
  shouldFetch: boolean
): UseQueryResult<DatasetUploadStatus> {
  const queryClient = useQueryClient();
  return useQuery<DatasetUploadStatus>(
    [USE_DATASET_STATUS, dataset_id],
    () => fetchDatasetStatus(dataset_id),
    {
      enabled: shouldFetch,
      onSuccess: async (data: DatasetUploadStatus): Promise<void> => {
        // When the dataset has been successfully processed, invalidate the collection query, and the datasets index query.
        // The collection query will fetch with the updated dataset list, which will no longer be in a loading state.
        // As a result, the useDatasetStatus query's status will be updated to "idle".
        if (data.processing_status === PROCESSING_STATUS.SUCCESS) {
          await queryClient.invalidateQueries([USE_COLLECTION, collectionId]);
          await queryClient.invalidateQueries(
            [USE_DATASETS_INDEX],
            DEFAULT_BACKGROUND_REFETCH
          );
        }
      },
      refetchInterval: REFETCH_INTERVAL_MS,
    }
  );
}

export const USE_DELETE_DATASET = {
  entities: [ENTITIES.DATASET],
  id: "dataset",
};

export interface DeleteDataset {
  collectionId: Collection["id"];
  datasetId: Dataset["id"];
}

async function deleteDataset({
  collectionId,
  datasetId,
}: DeleteDataset): Promise<DeleteDataset> {
  const url = apiTemplateToUrl(API_URL + API.DATASET, {
    dataset_id: datasetId,
  });
  const response = await fetch(url, DELETE_FETCH_OPTIONS);

  if (!response.ok) {
    throw Error(response.statusText);
  }

  return { collectionId, datasetId };
}

export function useDeleteDataset(): UseMutationResult<
  DeleteDataset,
  unknown,
  DeleteDataset
> {
  const queryClient = useQueryClient();
  return useMutation(deleteDataset, {
    onSuccess: async ({
      collectionId,
      datasetId,
    }: DeleteDataset): Promise<void> => {
      // If the dataset was in the process of loading, the dataset status query will be cancelled.
      // This is not an essential step, but cancels the request prior to invalidating the collection query (which
      // would have cancelled the dataset status query anyway).
      await queryClient.cancelQueries([USE_DATASET_STATUS, datasetId]);
      // Invalidate the collection query, and datasets index query.
      // Invalidation of the collection query triggers an immediate re-fetch with the dataset removed from the
      // list of datasets. Note, this re-fetch happens before the "useDeleteDataset" mutate function executes the
      // "onSuccess" callback.
      await queryClient.invalidateQueries([USE_COLLECTION, collectionId]);
      await queryClient.invalidateQueries(
        [USE_DATASETS_INDEX],
        DEFAULT_BACKGROUND_REFETCH
      );
    },
  });
}

/**
 * Query key for /datasets/id/assets.
 */
export const USE_DATASETS_ASSETS = {
  entities: [ENTITIES.DATASET],
  id: "datasetAsset",
};

/**
 * Cache-enabled hook for fetching assets for the dataset with the given ID.
 * @param datasetId - ID of dataset to fetch assets of.
 * @param enabled - True if fetch can be invoked.
 * @returns Dataset metadata.
 */
export function useFetchDatasetAssets(
  datasetId: string,
  enabled: boolean
): UseQueryResult<DatasetAsset[]> {
  return useQuery<DatasetAsset[]>(
    [USE_DATASETS_ASSETS, datasetId],
    () => fetchDatasetAssets(datasetId),
    { enabled }
  );
}

/**
 * Fetch assets for the dataset with the given ID.
 * @param datasetId - ID of dataset to fetch assets of.
 * @returns Promise that resolves to dataset metadata.
 */
async function fetchDatasetAssets(datasetId: string): Promise<DatasetAsset[]> {
  const { assets } = await (
    await fetch(
      apiTemplateToUrl(API_URL + API.DATASET_ASSETS, {
        dataset_id: datasetId,
      }),
      DEFAULT_FETCH_OPTIONS
    )
  ).json();
  return assets;
}
