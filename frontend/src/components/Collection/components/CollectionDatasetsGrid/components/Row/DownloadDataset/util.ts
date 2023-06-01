import { DatasetAsset } from "src/common/entities";
import { useFetchDatasetAssets } from "src/common/queries/datasets";

export interface UseDatasetAssets {
  datasetAssets: DatasetAsset[];
  isError: boolean;
  isLoading: boolean;
}

/**
 * Fetch dataset assets.
 * @param datasetId - ID of dataset to fetch dataset assets for.
 * @param enabled - True if fetch is enabled.
 * @returns Fetched dataset assets, if any, as well as fetch loading and error states.
 *
 */
export function useDatasetAssets(
  datasetId: string,
  enabled: boolean
): UseDatasetAssets {
  // Fetch dataset metadata on open of download modal.
  const {
    data: datasetAssets,
    isLoading,
    isError,
    isSuccess,
  } = useFetchDatasetAssets(datasetId, enabled);

  // Check for error state where dataset assets are missing - fetch is enabled, data is loaded but there are no
  // dataset assets.
  const isMissingAssets = isSuccess && !datasetAssets;

  return {
    datasetAssets: datasetAssets ?? [],
    isError: isError || isMissingAssets,
    isLoading,
  };
}
