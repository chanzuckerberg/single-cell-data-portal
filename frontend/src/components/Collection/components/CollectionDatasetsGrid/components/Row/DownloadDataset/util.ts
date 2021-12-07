import { Dataset, DatasetAsset } from "src/common/entities";
import { useFetchDatasetMetadata } from "src/common/queries/datasets";

export interface UseDatasetAssets {
  datasetAssets: DatasetAsset[];
  isError: boolean;
  isLoading: boolean;
}

/**
 * Fetch dataset assets.
 * @param datasetId - ID of dataset to fetch dataset assets for.
 * @param explorerUrl - Dataset's explorer URL, used to fetch dataset metadata.
 * @param enabled - True if fetch is enabled.
 * @returns Fetched dataset assets, if any, as well as fetch loading and error states.
 *
 */
export function useDatasetAssets(
  datasetId: string,
  explorerUrl: string,
  enabled: boolean
): UseDatasetAssets {
  // Fetch dataset metadata on open of download modal.
  const {
    data: datasetMetadata,
    isLoading,
    isError,
    isSuccess,
  } = useFetchDatasetMetadata(explorerUrl, enabled);

  // Find the dataset assets for the current dataset.
  const dataset = datasetMetadata?.collection_datasets.find(
    (collectionDataset: Dataset) => collectionDataset.id === datasetId
  );
  const datasetAssets = dataset?.dataset_assets;

  // Check for error state where dataset assets are missing - fetch is enabled, data is loaded but there are no
  // dataset assets.
  const isMissingAssets = isSuccess && !datasetAssets;

  return {
    datasetAssets: datasetAssets ?? [],
    isError: isError || isMissingAssets,
    isLoading,
  };
}
